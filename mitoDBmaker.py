# MitoDBmaker.py

'''	
main()
    Reads in each line of a file containing id,family,genus,species,subspecies with no spaces.
    It is recommended the list be aready ordered by family, genus, then species (checks for redundancy).
    If the mtGenome is present, gets all genes from one XML file,
    otherwise searches for separate gene XML files
get_ids()
    Gets giID (or list of IDs) from the nucleotide database
    This id is required for retrieval of the xml file in the getXML function
report_synonym()
    When the genus and species are listed under another name in the NCBI Nucleotide database,
    this function takes the genus and species names and enters it into a table as a synonym
    and uses the names it is "translated" to (transgenus and trans_species) as the reference.
get_full_taxonomy(genus,species)
    Given the Genus and Species of an organism, this function will search for the full taxonomy.
    This taxonomy will be entered as attributes (columns) in the Sqlite database.
getXML()
    Searches for and returns the XML file for the mitochondrial gene or genome from the nucleotide database.
    Checks that organism name is correct for genus and species.
get_gene_from_xml()
    Checks gene name, gets "to" and "from", takes compliment of sequence if neccessary, then extracts sequence
    Checks if the Genus, Species, Gene, and accession number and sequence already exist in the database
    updates if yes, inserts new row if no
retryResponse()
    If a bad conection occurs, this function will retry the connection a set number of times (e.g., 100)
    resting 0.3 seconds between each retry
Note: Currently ignores family that is given because I'm not sure if these are all the correct families
'''

from bs4 import BeautifulSoup  # to read XML
from ast import literal_eval  # for reading gene names from the terminal
from difflib import SequenceMatcher  # for finding the best match to a genes name (e.g. nad2 for ND2)
import sqlite3  # for search, insert, and update to database
import httplib  # to access NCBI URLs
import urllib
import string  # for string.replace in filename
import time
import sys
# import re  # for regular expressions

lastRequestTime = 0


def get_ids(genus, species, subspecies, gene_or_genome):
    global lastRequestTime
    #############################################################
    # get the GenBank accession number via esearch, ex 256557273
    #############################################################
    # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=Atractaspis+bibronii+COX1
    # Note: check that if subspecies isn't found, still finds species
    if subspecies is not None:
        url = '/entrez/eutils/esearch.fcgi?db=nuccore&term={}+{}+{}+{}'.format(
            urllib.quote(genus), urllib.quote(species), urllib.quote(subspecies),  
            urllib.quote(gene_or_genome))
    else:
        url = '/entrez/eutils/esearch.fcgi?db=nuccore&term={}+{}+{}'.format(
            urllib.quote(genus), urllib.quote(species), urllib.quote(gene_or_genome))

    #print("searching for " + gene_or_genome + " giID " + url)
    #NCBI requires a .3 second delay between requests
    time_since_last_request = time.time()-lastRequestTime
    if time_since_last_request<0.3:
        time.sleep(0.3-time_since_last_request)
    response = get_response(url)
    lastRequestTime = time.time()

    if response is None:
        raise ValueError('No xml data returned for {} {} {}. Query:\n{}'.format(genus, species, gene_or_genome, url))
    r_data = response.read()
    gi_soup = BeautifulSoup(r_data, 'xml')

    gi_id_list = gi_soup.IdList

    if gi_id_list is None:  # Genus not found, same gene_or_genome returned
        print "########## No id list found for " + str(genus) + " " + str(species) + " " + \
              str(gene_or_genome)
        return None, None, None

    translation_set = gi_soup.TranslationSet.Translation

    if translation_set is None:  # Genus not found, same gene_or_genome returned
        print "########## No result returned for " + str(genus) + " " + str(species) + " " + str(gene_or_genome)
        return None, None, None

    translation_string = translation_set.To.string
    translation = translation_string.split('"')
    translation = translation[1].split(' ')

    if len(translation) == 1:  # does not return genus and species names, only genus
        #print(" species " + species + " not available")

        return None, None, None

    else:
        trans_genus = translation[0]
        trans_species = translation[1]

    if gi_id_list.Id is None:  # No id returned for gene_or_genome but name translation available
        print "########## No giID for " + str(genus) + " " + str(species) + " " + \
              str(gene_or_genome)
        #print(" trans_genus is " + trans_genus + " and trans_species is " + trans_species)
        return None, trans_genus, trans_species

    # giIDs = gi_soup.IdList.Id.string #only returns first id
    gi_id_list = gi_soup.IdList
    giIDs = gi_id_list.find_all("Id")  # return a list of IDs

    return giIDs, trans_genus, trans_species  # returns "translated" name also


def report_synonym(genus, species, trans_genus, trans_species, db_connection):
    con = db_connection
    cur = con.cursor()

    organism = str(genus + " " + species)
    trans_organism = str(trans_genus + " " + trans_species)
    try:
        select = 'SELECT * FROM Synonyms WHERE Synonym=?'
        cur.execute(select, (organism,))

    except sqlite3.OperationalError:
        sql_command = '''CREATE TABLE "Synonyms" (
        "Synonym" TEXT PRIMARY KEY ON CONFLICT ABORT UNIQUE ON CONFLICT ABORT NOT NULL ON CONFLICT ABORT,
        "ReferenceName" TEXT NOT NULL ON CONFLICT ABORT
        ); '''
        cur.execute(sql_command)
        select = 'SELECT * FROM Synonyms WHERE Synonym=?'
        cur.execute(select, (organism,))

    rowcount = len(cur.fetchall())
    if rowcount < 0:
        print "invalid rowcount returned from database"
    # if rowcount == 1, Entry exists, but do not want to update Synonyms without updating RefGenes first
    if rowcount == 0: # insert
        insert = 'INSERT INTO Synonyms (Synonym,ReferenceName) VALUES (?,?)'
        cur.execute(insert, (organism, trans_organism))
    con.commit()

def get_full_taxonomy(genus, species):
    global lastRequestTime
    AllOtherRank = []
    Superkingdom = None
    Kingdom = None
    Superphylum = None
    Phylum = None
    Subphylum = None
    Class = None
    Subclass=None
    Superorder = None
    Order = None
    Suborder = None
    Infraorder = None
    Parvorder = None
    Superfamily = None
    Family = None
    Subfamily = None
    Tribe=None
    #####################################################################
    # get id from taxonomy database for specified genus and species
    #####################################################################
    # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=Hypsiglena+slevini
    if species!="sp.":
        url = '/entrez/eutils/esearch.fcgi?db=taxonomy&term=' + genus + '+' + species
    else:
        url = '/entrez/eutils/esearch.fcgi?db=taxonomy&term=' + genus
    print("searching for taxonomy " + url)


    time_since_last_request = time.time()-lastRequestTime
    if time_since_last_request<0.3:
        time.sleep(0.3-time_since_last_request)
    response = get_response(url)
    lastRequestTime = time.time()

    r_data = response.read()
    id_soup = BeautifulSoup(r_data, 'xml')

    if id_soup is None:
        print "No taxonomy xml returned"
        return

    if id_soup.eSearchResult.IdList.Id is None:
        print "Notaxonomy id for: " + str(genus) + " " + str(species)
        return

    id = id_soup.eSearchResult.IdList.Id.string

    #####################################################
    # Get taxonomy xml for the id
    #####################################################
    # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=241197&retmode=xml
    url = '/entrez/eutils/efetch.fcgi?db=taxonomy&id=' + id + '&retmode=xml'
    print("Getting taxonomy xml " + str(url))

    time_since_last_request = time.time()-lastRequestTime
    if time_since_last_request < 0.3:
        time.sleep(0.3-time_since_last_request)
    response = get_response(url)
    lastRequestTime = time.time()

    r_data = response.read()
    fastaSoup = BeautifulSoup(r_data, 'xml')
    

    ###################################################
    # check that the genus and species names are correct
    ###################################################

    org_list = fastaSoup.find("ScientificName")  # start value in list format
    organism = org_list.contents[0]  # Genus species
    org_name = str(organism).split()

    # lineage = fastaSoup.find("LineageEx")
    taxon = fastaSoup.find_all("Taxon")

    if org_name[0] != genus:
            correct_organism = False
            synonyms = taxon[0].find_all("Synonym")
            synonyms += taxon[0].find_all("GenbankSynonym")

            for synonym in synonyms:
                if synonym is None or synonym.string is None:
                    continue
                if len(org_name)<2:  # missing species name, genus only
                    # check for synonym
                    if genus in synonym.string:
                        correct_organism = True
                        break
                else:
                    scientific_name = genus + " " + species
                    if scientific_name in synonym.string:
                        correct_organism = True
                        break
                    else:
                        print(scientific_name + " not found in " + synonym)
            if correct_organism == False:
                print "No synonym found for " + str(genus) + " " + str(species)
                return


    ###############################################
    # Get the lineage
    ###############################################

    for taxa in taxon:

        rank = taxa.find("Rank")
        if rank.string == "no rank":
            AllOtherRank += [str(taxa.ScientificName.string)]
        elif rank.string == "superkingdom":
            Superkingdom = taxa.ScientificName.string
        elif rank.string == "kingdom":
            Kingdom = taxa.ScientificName.string
        elif rank.string == "superphylum":
            Superphylum = taxa.ScientificName.string
        elif rank.string == "phylum":
            Phylum = taxa.ScientificName.string
        elif rank.string == "subphylum":
            Subphylum = taxa.ScientificName.string
        elif rank.string == "class":
            Class = taxa.ScientificName.string
        elif rank.string == "subclass":
            Subclass = taxa.ScientificName.string
        elif rank.string == "superorder":
            Superorder = taxa.ScientificName.string
        elif rank.string == "order":
            Order = taxa.ScientificName.string
        elif rank.string == "suborder":
            Suborder = taxa.ScientificName.string
        elif rank.string == "infraorder":
            Infraorder = taxa.ScientificName.string
        elif rank.string == "parvorder":
            Parvorder = taxa.ScientificName.string
        elif rank.string == "superfamily":
            Superfamily = taxa.ScientificName.string
        elif rank.string == "family":
            Family = taxa.ScientificName.string
        elif rank.string == "subfamily":
            Subfamily = taxa.ScientificName.string
        elif rank.string == "genus":
            continue
        elif rank.string == "tribe":
            Tribe = taxa.ScientificName.string
        elif rank.string == "species":
            continue
        # Subspecies is taken from the organism name in the xml, usually None is given
        else:
            print "unknown rank " + str(rank.string) + " " + str(taxa.ScientificName.string)
            AllOtherRank += [str(taxa.ScientificName.string)]

        # if rank.string=="no rank" and taxa.ScientificName.string=="Cetartiodactyla":
            # order = "Cetartiodactyla"

    return str(AllOtherRank), Superkingdom, Kingdom, Superphylum, Phylum, Subphylum, Class, Subclass, Superorder, Order, \
           Suborder, Infraorder, Parvorder, Superfamily, Family, Subfamily, Tribe


def getXML(genus, species, gene_or_genome, giID):
    global lastRequestTime
    subspecies = None

    ####################################################
    # get sequence for specified gene or genome
    ####################################################
    # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=1043616945&retmode=xml
    # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=664680423&retmode=xml

    print(" getting sequence from " + giID + " xml")

    url = '/entrez/eutils/efetch.fcgi?db=nuccore&id=' + giID + '&retmode=xml'

    time_since_last_request = time.time()-lastRequestTime
    if time_since_last_request<0.3:
        time.sleep(0.3-time_since_last_request)
    response = get_response(url)  # get_full_taxonomy
    lastRequestTime=time.time()

    r_data = response.read()

    fasta_soup = BeautifulSoup(r_data, 'xml')
    # NOTE: BeautifulSoup does not allow you to parse dashes in xml names

    # check that genus and species name are correct
    org_list = fasta_soup.find("GBSeq_organism") #start value in list format
    # print "full organism name is " + str(org_list)
    organism = org_list.contents[0]  # Genus species
    org_name = str(organism).split()

    if genus != org_name[0] or species != org_name[1]:
        ##### check first two words of organism in case there is a subspecies listed ####
        print '##### Incorrect organism ' + str(organism) + ' for ' + str(genus) + ' ' + str(species) + ' ' +\
              str(gene_or_genome)
        return None, None  # means it failed

    #if gene_or_genome == 'mitochondrion, complete genome':  # save sequence
    if gene_or_genome == 'chloroplast, complete genome':  # save sequence
        geneLoc = fasta_soup.find("GBSeq_locus")
        geneLocus = geneLoc.contents[0]
        geneLocus = geneLocus.upper()  # accession number can be written in upper- or lowercase
        if geneLocus.startswith('NC') is False:
            print("geneLocus " + geneLocus + " is not a chromosome from RefSeq, Does it start with EU?")
            return None, None

    if len(org_name)>2:
        print "############### subspecies is " + str(org_name[2]) + " for " + str(genus) + " " + str(species)
        subspecies = str(org_name[2])

    return fasta_soup, subspecies


def get_gene_from_xml(fasta_soup, gene, db_connection, subspecies=None):  # fasta_Soup is an object that contins the XML with the gene

    # check for synonyms
    # matching ratio is too high when comparing "NADH 1" and "NADH 2",
    # so must use exact match for longer words
    # Lists originally taken from http://www.genecards.org/ which has multiple sources including HUGO


    gene_synonym = None
    '''
    if gene =="COX1":
        gene_synonym = ["CO1", "COI", "COXI", "Cytochrome C Oxidase I",
        "Cytochrome C Oxidase Subunit I", "MTCO1",
        "Mitochondrially Encoded Cytochrome C Oxidase I",
        "Cytochrome C Oxidase Polypeptide I", "EC 1.9.3.1"]
    if gene == "CYTB":
        gene_synonym = ["COB", "Cytochrome B",  "cyt b",
        "Mitochondrially Encoded Cytochrome B", "Complex III Subunit 3",
        "Complex III Subunit III", "Cytochrome B-C1 Complex Subunit 3",
        "Ubiquinol-Cytochrome-C Reductase Complex Cytochrome B Subunit", "MTCYB"]
    if gene == "ND2":
        # must search for exact match here, seqMatcher will match "NADH" with "NADH2"
        gene_synonym = ["NADH2", "NADH Dehydrogenase 2", "NADH Dehydrogenase Subunit 2",
        "Complex I ND2 Subunit", "NADH-Ubiquinone Oxidoreductase Chain 2", "MTND2",
        "Mitochondrially Encoded NADH Dehydrogenase 2", "EC 1.6.5.3"]
    if gene =="12S":  # "12S rRNA gene" AND "12S ribosomal RNA" should be captured by 12S alone
        gene_synonym = ["12S RNA", "s-rRNA", "Mitochondrially Encoded 12S RNA", "MTRNR1", "RNR1"] #short
        # note: gene for 12S is MT-RNR1 in animals
    if gene =="16S":
        gene_synonym = ["l-rRNA", "Mitochondrially Encoded 16S RNA","MTRNR2", "RNR2",
        "Humanin Mitochondrial", "Formyl-Humanin", "Humanin", "HNM", "HN"]
    #if gene == "trnV":
    #    gene_synonym = ["tRNA-Val","TRNA Valine","Mitochondrially Encoded TRNA Valine", "MTTV"]
    #    # Note: TRNA is a known synonym but there are other types of TRNA. Check for this!!!
    if gene == "ND5":
        gene_synonym = ["Mitochondrially Encoded NADH:Ubiquinone Oxidoreductase Core Subunit 5",
        "NADH Dehydrogenase Subunit 5","EC 1.6.5.3", "MTND5", 
        "Mitochondrially Encoded NADH Dehydrogenase 5", 
        "NADH Dehydrogenase, Subunit 5 (Complex I)", "NADH-Ubiquinone Oxidoreductase Chain 5",
        "Complex I ND5 Subunit","NADH Dehydrogenase 5","NADH5"]
    '''
    if gene == 'rpoB':
        gene_synonym = ['RNA polymerase beta subunit']
    if gene == 'rpoC1':
        gene_synonym = ['RNA polymerase C']
    if gene == 'trnC':
        gene_synonym = ['trnC(gca)']
    if gene == 'trnG':
        gene_synonym = ['tRNA-Gly', 'tRNA-Gly (trnG)']
    if gene == 'trnE':
        gene_synonym = ['tRNA-Glu (trnE)']
    if gene == 'trnH':
        gene_synonym = ["tRNA-His"]
    if gene == 'trnK':
        gene_synonym = ['trnK-UUU', 'trnK(ctt)']
    if gene == 'trnL':
        gene_synonym = ['trnL(UAA)']
    if gene == 'trnS':
        gene_synonym = ['trnS-GCU', 'trnS-GGA']
    if gene == 'ITS1':
        gene_synonym = ['internal transcribed spacer 1', 'internal transcribed spacer 1, ITS1']
    if gene == 'ITS2':
        gene_synonym = ['internal transcribed spacer 2', 'internal transcribed spacer 2, ITS2']
    if gene == 'ITS':
        gene_synonym = ['sequence contains ITS1, 5.8S rRNA gene, ITS2', 'sequence contains 18S rRNA gene, ITS1, 5.8S rRNA gene, ITS2, 28S rRNA gene']

    if gene == 'atpB+rbcL':
        gene_synonym = ['atpB-rbcL', 'atpB-rbcL intergenic spacer']
    if gene == 'atpF+atpH':
        gene_synonym = ['atpF-atpH', 'atpF-atpH intergenic spacer']
    if gene == 'ndhF+rpl32':
        gene_synonym = ['ndhF-rpl32', 'ndhF-rpl32 intergenic spacer']
    if gene == 'ndhJ+trnC':
        gene_synonym = ['ndhJ-trnC','ndhJ-trnC intergenic spacer']
    if gene == 'psbA+trnH':
        gene_synonym = ['psbA-trnH', 'psbA-trnH intergenic spacer', 'psbA-trnH intergenic spacer region'
                        'contains portions of discontiguous regions including ribosomal protein L16 (rpl16) intron, psbA-trnH intergenic spacer, and trnC-ycf6 intergenic spacer'
                        'psbA-trnH intergenic spacer; may also contain tRNA-His (trnH)']
    if gene == 'psbD+trnT':
        gene_synonym = ['psbD-trnT', 'psbD-trnT intergenic spacer']
    if gene == 'psbK+psbI':
        gene_synonym = ['psbK-psbI', 'psbK-psbI intergenic spacer', 'contains psbK-psbI intergenic spacer']
    if gene== 'psbZ+trnG':
        gene_synonym = ['psbZ-trnG','psbZ-trnG intergenic spacer and tRNA-Gly',
                        'contains psbZ-trnG intergenic spacer and tRNA-Gly']

    if gene == 'rpl32+trnL':
        gene_synonym = ['rpl32-trnL', 'contains rpl32 gene (partial)', 'rpl32-trnL(UAG) intergenic spacer'
                        'rpl32-trnL IGS and trnL gene (partial)',
                        'contains rpl32 gene and rpl32-trnL intergenic spacer',
                        'contains rpl32, rpl32-trnL intergenic spacer and tRNA-Leu',
                        'rpl32-trnL intergenic spacer','rpl32-trnL intergenic spacer region']



    if gene == 'rp132+ndhF':
        gene_synonym = ['rp132-ndhF','rpl32-ndhF intergenic spacer']
    if gene == 'rpoB+trnC':
        gene_synonym = ['rpoB-trnC', 'rpoB-trnC intergenic spacer']
    if gene == 'rps4+trnS':
        gene_synonym = ['rps4-trnS', 'rps4-trnS intergenic spacer (trnS-GGA)']
    if gene == 'rps16+trnQ':
        gene_synonym = ['rps16-trnQ', 'rps16-trnQ intergenic spacer region']
    if gene == 'rps16+trnK':
        gene_synonym = ['rps16-trnK', 'rps16-trnK intergenic spacer']
    if gene == 'trnC+ycf6':
        gene_synonym = ['trnC-ycf6','trnC-ycf6 intergenic spacer', 'trnC-yCf6 intergenic spacer region',
                        ' tRNA-Cys (trnC) and trnC-ycf6', 'tRNA-Cys (trnC) and trnC-ycf6 intergenic spacer']
    if gene== 'trnC+petN':
        gene_synonym = ['trnC-petN', 'trnC-petN intergenic spacer, IGS', 'trnC-petN intergenic spacer']
    if gene == 'trnD+trnT':
        gene_synonym = ['trnD-trnT', 'trnD-trnT intergenic spacer']
    if gene == 'trnD+trnY':
        gene_synonym == ['trnD-trnY', 'contains tRNA-Asp (trnD), trnD-trnY intergenic spacer, tRNA-Tyr (trnY), trnY-trnE intergenic spacer, tRNA-Glu (trnE), trnE-trnT intergenic spacer, and tRNA-Thr (trnT)']

    if gene == 'trnE+trnT':
        gene_synonym = ['trnE-trnT','trnE-trnT intergenic spacer', 'chloroplast marker trnET; chloroplast trnE/trnT intergenic region',
                        'tRNA-Asp (trnD), trnD-trnY intergenic spacer, tRNA-Tyr (trnY), trnY-trnE intergenic spacer, tRNA-Glu (trnE), trnE-trnT intergenic spacer, and tRNA-Thr (trnT)']
    if gene == 'trnE+trnY':
        gene_synonym =  ['trnE-trnY', 'trnE-trnY intergenic spacer',
                         'trnD-trnY intergenic spacer, tRNA-Tyr (trnY), trnY-trnE intergenic spacer, tRNA-Glu (trnE), trnE-trnT intergenic spacer, and tRNA-Thr (trnT)']
    if gene == 'trnG+trnS':
        gene_synonym = ['trnG-trnS','trnG-trnG intergenic spacer',
                        'tRNA-Gly (trnG) gene and trnG-trnS intergenic spacer']

    if gene == 'trnH+psbA':
        gene_synonym = ['trnH-psbA', 'trnH-psbA intergenic spacer', 'trnH-psbA intergenic spacer region']
    if gene == 'trnK+matK':
        gene_synonym= ['trnK-matK intergenic spacer']

    if gene == 'trnL+trnF':
        gene_synonym = ['trnL-trnF','contains tRNA-Leu, trnL-trnF intergenic spacer and tRNA-Phe', 'tRNA-Leu', 'trnL-trnF intergenic spacer', 'contains tRNA-Leu (trnL) and trnL-trnF intergenic spacer'
                        'contains tRNA-Leu (trnL), trnL-trnF intergenic spacer, and tRNA-Phe (trnF)']
    if gene == 'trnQ+rps16':
        gene_synonym = ['trnQ-rps16','trnQ(UUG)-rps16 intergenic spacer']
    if gene == 'trnS+psbZ':
        gene_synonym =  ['trnS-psbZ', 'tRNA-Ser and trnS-psbZ intergenic spacer', 'contains tRNA-Ser and trnS-psbZ intergenic spacer']
    if gene == 'trnS+trnfM':
        gene_synonym = ['trnS-trnfM', 'trnS-trnfM intergenic spacer']
    if gene == 'trnS+trnG':
        gene_synonym= ['trnS-trnG', 'trnS-trnG intergenic spacer', 'trnS-trnG intergenic spacer region', 'tRNA-Ser (trnS), trnS-trnG intergenic spacer and tRNA-Gly (trnG)']
    if gene == 'trnT+trnL':
        gene_synonym = ['trnT-trnL', 'trnT-trnL', 'trnT-trnL intergenic spacer', 'trnT-trnL intergenic spacer region',
                        'trnT-trnF intergenic spacer region', 'trnT(UGU)-trnL(UAA) intergenic spacer', 'contains trnT-trnL intergenic spacer, trnL gene and trnL-trnF intergenic spacer sequences', 'trnT-trnL intergenic spacer, tRNA-Leu gene, and trnL-trnF intergenic spacer', 'contains trnT-trnL intergenic spacer, tRNA-Leu (trnL), trnL-trnF intergenic spacer, and tRNA-Phe (trnF)']
    if gene == 'trnY+trnE':
        gene_synonym = ['trnY-trnE', 'trnY-trnE intergenic spacer','trnD, trnY, trnE, and trnE-trnT intergenic spacer region']
    if gene == 'ycf6+psbM':
        gene_synonym = ['ycf6-psbM', 'ycf6-psbM intergenic spacer']
    if gene == 'ycf6+trnC':
        gene_synonym = ['ycf6-trnC','ycf6-trnC-GCA intergenic spacer']




    gene_loc = fasta_soup.find("GBSeq_locus")
    gene_locus = gene_loc.contents[0]
    # insert gene_locus into database, this is the accession number

    genus_species = fasta_soup.find("GBSeq_organism")
    genus_species = genus_species.contents[0]
    genus_species = genus_species.split()
    genus = genus_species[0]
    species = genus_species[1]
    if str(species) == 'x':
        print "hybrid found"
        return None

    try:
        subspecies = genus_species[2]
    except IndexError:
        subspecies = None

    seqList = fasta_soup.find("GBSeq_sequence")

    if seqList is None:  # File contains only GBSeq_contig
        print "		The xml file returned has no sequence. May be a file listing contigs."
        return None
    full_sequence = seqList.contents[0]  # full sequence

    # Find closest match to gene we are looking for
    exact_match_found = False
    best_match = 0  # Start with no match, search for best
    best_feature = None
    best_word = None

    features = fasta_soup.find_all("GBFeature")

    for feat in features:
        if exact_match_found == True:
            break
        if feat.GBFeature_quals is None:
            continue
        # qualifiers = feat.GBFeature_quals.GBQualifier only finds first
        qualifiers = feat.find_all("GBQualifier")
        for qualifier in qualifiers:
            if qualifier.GBQualifier_value is None:
                continue
            value = qualifier.GBQualifier_value.string
            #Sprint("value is " + value)
            # Check entire value for a match before individual words
            if value.upper()==gene.upper():
                print(" ##### exact match #### " + str(gene))
                exact_match_found=True
                best_match = 1.0
                best_feature = feat
                best_word = value
                break

            if gene_synonym is not None:  # look for exact match of synonym
                #print(value.upper() + " did not match " + gene.upper())
                for syn in gene_synonym:
                    if value.upper() == syn.upper():
                            best_match = 1.0
                            best_feature = feat
                            best_word = value
                            print "##### gene synonym match #### " + str(syn)
                            break  # want to break out of outer loop, make function?

            # for each word in value
            # check the matching ratio
            words = value.split()
            # print "words are " + str(words)

            # if no exact match, check individual words
            for word in words:
                s = SequenceMatcher(None, word.upper(), gene.upper())
                matchRatio=s.ratio()  # 0 is no match, 1 is perfect match
                # if ratio is better that best match, save as best_feature and make new best_match score
                if matchRatio>best_match:
                    best_match=matchRatio
                    best_feature=feat
                    best_word=word

    if best_feature is None:
        print str(gene) + " not found in XML file for " + genus + " " + species
        return None

    if best_match is None:
        print str(gene) + " has no match found in XML file for " + genus + " " + species
        return None

    gene_match = best_feature.GBFeature_quals.GBQualifier.GBQualifier_value.string

    if 1 > best_match > 0.8:  # else below (prints closest) or perfect match (no print)
        print "closest match for " + gene + " is " + str(gene_match)

    if best_match <= 0.8:
        # print "############### No close match: " + str(best_word)
        # print " for " + str(gene) + " of " + str(genus) + ' ' + str(species)
        print "############### No close match: " + str(gene_match) + " for " + str(gene) + " of " + str(genus) + ' ' \
              + str(species) + ' ' + str(gene_locus)

        return None

    # get start and stop
    location = best_feature.GBFeature_location.string
    #has_join = False
    sequence=""
    #Check for a join statement

    if "join" in location:
        # check for complement()around the entire join and remove
        is_complement = location.startswith("complement")
        #remove complement, join, and all parentheses
        if is_complement:
            location = location.strip("complement()")
        # check for join()around the entire location and remove
        location = location.strip("join()")
        locations = location.split(",")
        # retrieve location(start and stop) for each section and determine the number of nucleotides between (intron)
        i=0
        old_stop=0
        #print("locations are " + str(locations))
        for loc in locations:
            start, stop = loc.split("..")
            start, stop, exon, partial = get_sequence(start, stop, full_sequence, is_complement)
            sequence += exon

            if i==0: #if first in sequence, no intron
                intron = ""
            else:
                intron_length = int(start) - int(old_stop)
                intron = "x"*intron_length
            sequence += intron
            old_stop = stop
            # find number of joins by counting the number of commas
            i += 1

    else:
        start, stop = location.split("..")
        start, stop, sequence, partial = get_sequence(start, stop, full_sequence)

    # Save whether it is a gene or cds, rna, etc
    feature_type = best_feature.GBFeature_key.string

    # save organism for header, and insert as a field or column in database
    # organism = fasta_soup.GBSeq_organism.string
    '''
    fastaHeader = organism + "\t" + gene + "\t" + feature_type + "\t" + "partial=" + str(partial)
    #save fasta file
    #filePrefix = gene + '_' + organism.replace(" ","_")
    #f_out= open(filePrefix + ".fasta", 'w')
    #fasta_lines=(fastaHeader, '\n', sequence)
    #f_out.writelines(fasta_lines)
    #f_out.close()
    '''




    # submitToDatabase:
    # genus=>Genus, species=>Species, gene=>GeneName , gene_locus=>AccessionNumber, partial=>PartialFlag,
    # and sequence=>GeneSequence
    con = db_connection
    cur = con.cursor()

    # Check if the Genus, Species, Gene, and accession number and sequence already exist in the database
    # If it returns something, print the thing and do not update
    # otherwise update

    select = 'SELECT * FROM RefGenes WHERE Genus=? AND Species=? AND GeneName=? AND AccessionNumber=?'
    try:
        cur.execute(select, (genus, species, gene, gene_locus))
    # figure out how many rows it returned.

    except sqlite3.OperationalError, e:
        print str(e)
        print("Database table not found, creating tables.")
        sql_command = '''CREATE TABLE "RefGenes" (
        "Row" INTEGER NOT NULL ON CONFLICT ABORT PRIMARY KEY AUTOINCREMENT UNIQUE,
        "AllOtherRank" TEXT,
        "Superkingdom" TEXT,
        "Kingdom" TEXT,
        "Superphylum" TEXT,
        "Phylum" TEXT,
        "Subphylum" TEXT,
        "Class" TEXT,
        "Subclass" TEXT,
        "Superorder" TEXT,
        "Order" TEXT,
        "Suborder" TEXT,
        "Infraorder" TEXT,
        "Parvorder" TEXT,
        "Superfamily" TEXT,
        "Family" TEXT,
        "Subfamily" TEXT,
        "Genus" TEXT NOT NULL ON CONFLICT ABORT,
        "Tribe" TEXT,
        "Species" TEXT NOT NULL ON CONFLICT ABORT,
        "Subspecies" TEXT,
        "GeneName" TEXT NOT NULL ON CONFLICT ABORT,
        "GeneType" TEXT NOT NULL ON CONFLICT ABORT,
        "AccessionNumber" TEXT NOT NULL ON CONFLICT ABORT,
        "PartialFlag" INTEGER NOT NULL ON CONFLICT ABORT,
        "GeneSequence" BLOB NOT NULL ON CONFLICT ABORT,
        CONSTRAINT "NonDuplication" UNIQUE ("Genus", "Species", "GeneName", "AccessionNumber") ON CONFLICT ABORT
        );'''
        cur.execute(sql_command)
        
        sql_command = '''CREATE TABLE "Synonyms" (
        "Synonym" TEXT PRIMARY KEY ON CONFLICT ABORT UNIQUE ON CONFLICT ABORT NOT NULL ON CONFLICT ABORT,
        "ReferenceName" TEXT NOT NULL ON CONFLICT ABORT
        ); '''
        cur.execute(sql_command)


    rowcount = len(cur.fetchall())

    if rowcount < 0:
        print "invalid rowcount"
        return None

    if rowcount>0:  # update
        print "updating table for " + str(genus_species[0]) + ' ' + str(genus_species[1]) + ' ' + str(gene)
        #print get_full_taxonomy(genus, species)
        AllOtherRank, Superkingdom, Kingdom, Superphylum, Phylum, Subphylum, Class, Subclass, Superorder, Order, Suborder, \
        Infraorder, Parvorder, Superfamily, Family, Subfamily, Tribe = get_full_taxonomy(genus, species)

        AllOtherRank = AllOtherRank.replace("\'", "")
        update = ''' UPDATE RefGenes SET
            AllOtherRank=?, Superkingdom=?, Kingdom=?, Superphylum=?, Phylum=?, Subphylum=?, Class=?, Subclass=?,
            Superorder=?, "Order"=?, Suborder=?, Infraorder=?, Parvorder=?,
            Superfamily=?, Family=?, Subfamily=?, Genus=? ,  Tribe=?, Species=?, subspecies=?,
            GeneName=? , GeneType=?, AccessionNumber=? , PartialFlag=?, GeneSequence=?
            WHERE Genus=? AND Species=? AND GeneName=? AND AccessionNumber=? AND GeneSequence=?
            '''

        cur.execute(update, (AllOtherRank, Superkingdom, Kingdom, Superphylum, Phylum, Subphylum, Class, Subclass,
                             Superorder, Order, Suborder, Infraorder, Parvorder, Superfamily, Family, Subfamily, genus,
                             Tribe, species, subspecies, gene, feature_type, gene_locus, int(partial), sequence, genus,
                             species, gene, gene_locus, sequence))

    if rowcount == 0:  # insert
        try:
            AllOtherRank, Superkingdom, Kingdom, Superphylum, Phylum, Subphylum, Class, Subclass, Superorder, Order, Suborder, \
            Infraorder, Parvorder, Superfamily, Family, Subfamily, Tribe = get_full_taxonomy(genus, species)


        except:
            return None

        AllOtherRank = AllOtherRank.replace("\'", "")


        insert = 'INSERT INTO RefGenes (AllOtherRank, Superkingdom, Kingdom, Superphylum, Phylum, Subphylum, Class, ' \
                 'Subclass, Superorder, "Order", Suborder, Infraorder, Parvorder, Superfamily, Family, Subfamily,  ' \
                 'Genus, Tribe, Species, subspecies, GeneName, GeneType, AccessionNumber, PartialFLag, GeneSequence) ' \
                 'VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'
        cur.execute(insert, (AllOtherRank, Superkingdom, Kingdom, Superphylum, Phylum, Subphylum, Class, Subclass,
                             Superorder, Order,Suborder,Infraorder,Parvorder,Superfamily,Family,Subfamily, genus, Tribe,
                             species, subspecies, gene, feature_type, gene_locus, int(partial), sequence))

    con.commit()

    return "successful"

def get_sequence(start, stop, full_sequence, is_complement = False):

    partial = False
    if "complement" in start:  # check for complement() and remove
        start = start.strip("complement(")
        stop = stop.strip(")")
        is_complement = True

    if "<" in start or ">" in stop:  # check for &lt, &gt which indicates partial, and remove
        partial = True
        # insert partial into database
        start = start.strip("<")
        stop = stop.strip(">")

    start = int(start) - 1  # NCBI uses 1 for first item, python requires zero
    stop = int(stop)
    sequence = full_sequence[start:stop]

    # Complements require a==>t, g==> c, etc.
    translation_table = string.maketrans("atugcyrswkmbdhvn", "taacgryswmkvhdbn")

    if is_complement:
        translated_sequence = string.translate(str(sequence), translation_table)  # take complement
        sequence = translated_sequence[::-1]  # reverse order
    return start, stop, sequence, partial


def get_response(url):
    _server = 'eutils.ncbi.nlm.nih.gov'
    conn = httplib.HTTPSConnection(_server) 	# Setup HTTP session with the server
    conn.request('GET', url)  # get a connection
    response = conn.getresponse()

    if response.status == 200:
        return response

    else:
        retries = 100
        for i in range(retries):
            try:
                conn = httplib.HTTPSConnection(_server)
                conn.request('GET', url)
                response = conn.getresponse()
                if response.status == 200:
                    break
                time.sleep(.3)
            except (httplib.HTTPException):
                print "HTTPException in httplib"
        if response.status != 200:
            print('esearch Bad HTTP response ({}), retried 10 times and failed'.format(
                response.status))
            return None

    return response

if __name__ == '__main__':

    #######################################
    #read in list of Genus species names
    #######################################
    # Get the genus and species name from comma delimited file (e.x. ZaherRaw.csv)

    try:
        sample_list = sys.argv[1]  # raw data list from Alan's lab [id, family, genus, species]

    except IndexError:
        print "Using Default Filename"
        sample_list = "SampleList.csv"

    try:
        geneList = literal_eval(sys.argv[2])  # i.e. [\'COX1\', \'ND2\',\'12S\',\'16S\',\'ND5\']
        #Check that genes are qualified geneList names, to avoid redundant DB entries
        for gene in geneList:
            if gene.upper() == 'TRNG':
                print(' If trnG is not available. Try trnS as a seed.')
            if gene.upper() == "TRNL":
                print('If trnL not available. Try rpl32 or trnT as a seed.')
            if gene.upper() ==  'YCF6':
                print('If ycf6 is not available, try trnC as a seed.')



    except IndexError:
        print "Using default gene list"
        # VERTEBRATE GENE LIST
        #geneList = ['COX1', 'CYTB', 'ND2', '12S', '16S', 'ND5']

        #Angiosperm gene list
        # trnP(GGG)-psaJ
        '''
        geneList = ['accD', 'atpA', 'atpB', 'atpE', 'atpF', 'cemA', 'clpP', 'ccsA', 'infA', 'matK',
                    'ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhE', 'ndhG', 'ndhH', 'ndhI', 'ndhF', 'ndhJ', 'ndhK',
                    'petA', 'petB', 'petD',
                    'psaA', 'psaB', 'psaC', 'psaJ', 'psbA',
                    'psbB', 'psbC', 'psbD', 'psbE', 'psbF', 'psbH', 'psbK', 'psnL', 'psbN', 'psbT', 'psbZ',
                    'rbcL', 'rpl2', 'rpl12', 'rpl14', 'rpl16', 'rpl20', 'rpl23', 'rpl33',
                    'rpoC1', 'rpoC2', 'rpoA', 'rpoB', 'rpl2',
                    'rps2', 'rps4', 'rps7', 'rps8', 'rps11', 'rps12', 'rps14', 'rps15', 'rps16', 'rps18', 'rps19',
                    'rrn16', 'rrn23',
                    'trnA', 'trnC', 'trnE', 'trnG', 'trnH', 'trnH+psbA', 'trnI','trnK', 'trnL', 'trnS', 'trnT', 'trnY',
                    'ycf1', 'ycf2', 'ycf3', 'ycf4', 'ycf6',
                    'ITS1', 'ITS2', 'atpB+rbcL',  'atpF+atpH', 'ndhF+rpl32','ndhJ+trnC',
                    'psbA+trnH', 'psbD+trnT','psbK+psbI',  'psbZ+trnG', 'rpl32+trnL', 'rp132+ndhF',
                    'rpoB+trnC', 'rps4+trnS', 'rps16+trnQ', 'rps16+trnK', 'trnC+ycf6', 'trnC+petN',  'trnD+trnT',
                    'trnE+trnT', 'trnE+trnY', 'trnG+trnS', 'trnL+trnF', 'trnQ+rps16', 'trnS+psbZ',
                    'trnS+trnfM', 'trnS+trnG', 'trnT+trnL', 'trnY+trnE',  'ycf6+psbM', 'ycf6+trnC']
        '''
        #Plant Gene List
        geneList = ['atpB', 'atpB+rbcL', 'atpE', 'matK', 'ndhF', 'psbA+trnH',
                    'psbB', 'psbH', 'psbN', 'psbZ', 'psbZ+trnG',
                    'rbcL', 'rpl16', 'rpl32+trnL', 'rpoC1', 'rps4', 'rps16', 'trnD+trnY', 'trnK+matK',
                    'trnL+trnF', 'trnT+trnL', 'trnS+psbZ', 'ycf1', 'ITS1', 'ITS2']
        
        #Barcode gene list for angiosperms
        #geneList = ['rbcl', 'matK', 'ITS1', 'ITS2', 'ycf1', 'trnH+psbA', 'psbA+trnH']

        '''
        #Angiosperm gene List
        geneList = ['accD', 'atpA', 'atpB', 'atpE', 'atpF', 'cemA', 'clpP', 'ccsA', 'infA',
                    'ndhA', 'ndhB', 'ndhC', 'ndhD', 'ndhE', 'ndhG', 'ndhH', 'ndhI', 'ndhF', 'ndhJ', 'ndhK',
                    'matK', 'petA', 'petB', 'petD',
                    'psaA', 'psaB', 'psaC', 'psaJ',
                    'psbA', 'psbA+trnH',
                    'psbB', 'psbC', 'psbD', 'psbE', 'psbF', 'psbH', 'psbH', 'psbK', 'psnL', 'psbN', 'psbT', 'psbZ',
                    'rbcL', 'rbcLa', 'rpl2', 'rpl12', 'rpl14', 'rpl16', 'rpl20', 'rpl23', 'rpl33',
                    'rpoC1', 'rpoC2', 'rpoA', 'rpoB', 'rpl2',
                    'rps2', 'rps4', 'rps7', 'rps8', 'rps11', 'rps12', 'rps14', 'rps15', 'rps16', 'rps18', 'rps19',
                    'rrn16', 'rrn23',
                    'trnA', 'trnC', 'trnE', 'trnG', 'trnH', 'trnH+psbA', 'trnI','trnK', 'trnL', 'trnS', 'trnT', 'trnY',
                    'ycf1', 'ycf2', 'ycf3', 'ycf4', 'ycf6', 'ITS1', 'ITS2']
        '''


        '''
        #Intergenic spacers for angiosperms
        geneList = ['atpB+rbcL',  'atpF+atpH', 'ndhF+rpl32', 'ndhJ+trnC',  'psbA+trnH', 'psbD+trnT',
                    'psbK+psbI',  'psbZ+trnG', 'rpl32+trnL', 'rp132+ndhF', 'rpoB+trnC', 'rps4+trnS', 'rps16+trnQ',
                    'rps16+trnK', 'trnC+ycf6', 'trnC+petN',  'trnD+trnT',  'trnE+trnT', 'trnE+trnY', 'trnG+trnS',
                    'trnH+psbA',  'trnL+trnF', 'trnQ+rps16', 'trnS+psbZ',  'trnS+trnfM', 'trnS+trnG', 'trnT+trnL',
                    'trnY+trnE',  'ycf6+psbM', 'ycf6+trnC']
        '''

    old_genus = None
    old_species = None

    with open(sample_list, 'r') as sample:
        lines = sample.readlines()

        for line in lines:
            line = line.replace('\n', '')
            cells = line.split(',')
            ID = cells[0].strip()
            family = cells[1].strip()
            genus = cells[2].strip()
            species = cells[3].strip()

            try:
                subspecies = cells[4]
            except(IndexError):
                subspecies = None

            if subspecies == '' or subspecies == ' ':  # An extra comma on the end creates a subspecies
                subspecies = None

            outFileName = string.replace(sample_list, ".csv", "NeedReference.csv")
            fout = open(outFileName,'a')  # a for append

            if species == "" or species == "sp.":
                print "No species, writing to file"
                for gene in geneList:
                    fout.write(str(gene) + ',' + str(line) + "\n")
                continue

            # Redundancy check: the genus, in case same species as the previous line
            if genus == old_genus and species == old_species:
                # print "redundancy"
                continue
                
            dbConnection = sqlite3.connect('./mitoDB.sqlite')

            # Check for tables in database

            # If tables aren't there, make them
            # RefGenes
            # RefGenomes
            # Synonyms

            #genome='mitochondrion, complete genome'
            genome='chloroplast, complete genome'

            giIDs, trans_genus, trans_species = get_ids(genus, species, subspecies, genome)
            # get_ids checks for translated names (synonyms)
            #print "transgenus and transspecies are " + str(trans_genus) + " " + str(trans_species)

            if trans_genus is not None and trans_species is not None:
                if genus != trans_genus or species != trans_species:
                    # "translates" to another name in the database, original name becomes a synonym
                    report_synonym(genus, species, trans_genus, trans_species, dbConnection)
                    # synonym,synonym,reference,reference
                    genus = trans_genus  # set genus and species to the names given in the xml
                    species = trans_species

            if giIDs is not None:
                print("gi ids FOUND")
                giID = giIDs[0].string
                # Subspecies is whatever is listed in the database, may be different from the subspecies given in the sample list
                xml, Subspecies = getXML(genus, species, genome, giID)  # try getting XML for complete genome

            else:  # No id was returned for the genome, can't retrieve xml
                #print("gi id not found for complete genome ")
                xml = None

            entry = None  # "successful" if a sequence is entered into the Database

            if xml is not None:  # Genome present
                print(" xml found for complete genome ")
                for gene in geneList:  # get all genes from one XML file

                    entry = get_gene_from_xml(xml, gene, dbConnection, Subspecies)

                    if entry is None:
                        print("# gene doesn't match or gene location not given in xml")
                        # shouldn't happen unless missing a synonym
                        print str(gene) + " not found in genome of " + str(genus) + ' ' + str(species)
                        fout.write(str(gene) + ',' + str(line))
                    else:
                        print(gene + " gene found")


            else: # Genome not present, xml was None for Genome search, Now search for genes
                for gene in geneList:
                    giIDs, trans_genus, trans_species = get_ids(genus, species, subspecies, gene)

                    # trans_genus and trans_species already set during genome search, ignore them here
                    gene_entry = None

                    if giIDs is not None:
                        for giID in giIDs:
                            giID = giID.string
                            xml, Subspecies = getXML(genus, species, gene, giID)

                            if xml is not None:
                                gene_entry = get_gene_from_xml(xml, gene, dbConnection, Subspecies)

                                if gene_entry is not None:
                                    break 	# successful entry, no need to check the next id

                    if gene_entry is None:	# Most will end here, will search for genus-only in another program
                        fout.write(str(gene) + ',' + str(line) + "\n")
            old_genus = genus
            old_species = species

            fout.close()

            # Write function that populates relatives

            dbConnection.close()

    # We can determine closely related genus or family for the ones printed to file
    # and insert it if one does not already exist
    # closestMatches can come from mitoRefmatcher.csv, but this only works for whole mtGenomes

    #Ramphotyphlops,braminus error: 23 bindings supplied?