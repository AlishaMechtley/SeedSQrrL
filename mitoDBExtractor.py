# mitoDBExtractor.py
# Reads in each line of a file containing id,family,genus,species,subspecies with no spaces.
# For a each organism, finds the closest relative in the database
# Check if genus exists in mitoDB.sqlite, Do same for Subfamily, etc.
# Does not look in NCBI database for closest relative, will make another program to do this.

################ Note ############################
# consider possibly saving the full taxonomy for an organism
# so we don't have to look it up twice
# and maybe save the closest reference in a new table as well?
# or will this be a lot of extra info we don't want in the end


from collections import OrderedDict
from bs4 import BeautifulSoup
from ast import literal_eval
import sqlite3
import httplib
from mitoDBmaker import get_response
import sys

dbConnection = sqlite3.connect('mitoDB.sqlite')
dbConnection.row_factory = sqlite3.Row


def checkSynonyms(genus, species):
    trans_genus = None
    trans_species = None
    cur = dbConnection.cursor()
    organism = str(genus + " " + species)
    try:
        select = 'SELECT * FROM Synonyms WHERE Synonym=?'
        cur.execute(select, (organism,))
        fetched = cur.fetchone()
    except sqlite3.OperationalError:
        sql_command = '''CREATE TABLE "Synonyms" (
        "Synonym" TEXT PRIMARY KEY ON CONFLICT ABORT UNIQUE ON CONFLICT ABORT NOT NULL ON CONFLICT ABORT,
        "ReferenceName" TEXT NOT NULL ON CONFLICT ABORT
        ); '''
        cur.execute(sql_command)
        select = 'SELECT * FROM Synonyms WHERE Synonym=?'
        cur.execute(select, (organism,))
        fetched = cur.fetchone()
    if fetched is not None:  # not a synonym
        genus, species = str(fetched["Synonym"]).split()
        trans_genus, trans_species = str(fetched["ReferenceName"]).split()
        print(str(genus) + " translated to " + str(trans_genus) + " and " + str(species) + " translated to " + str(
            trans_species))

        # trans_genus=trans_genus.replace("u'", "")

    return trans_genus, trans_species


def get_taxonomy(rank_value, previous_rank_value):  # different return than DBMaker, uses dict; returns 3 fewer results
    print("getting taxonomy for " +  rank_value + " " + previous_rank_value)
    #print "rank_value is " + str(rank_value)
    #print "previous_rank_value is " + str(previous_rank_value)


    # AllOtherRank=[]
    Superkingdom = None
    Kingdom = None
    Superphylum = None
    Phylum = None
    Subphylum = None
    Class = None
    Subclass = None
    Superorder = None
    Order = None
    Suborder = None
    Infraorder = None
    Parvorder = None
    Superfamily = None
    Family = None
    Subfamily = None
    Tribe = None

    #####################################################################
    # get id from taxonomy database for specified genus and species
    #####################################################################
    # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=taxonomy&term=Eumetopias+jubatus
    url = '/entrez/eutils/esearch.fcgi?db=taxonomy&term=' + str(rank_value) + '+' + str(previous_rank_value)

   #print("taxonomy url is " + str(url))
    response = get_response(url)
    rData = response.read()
    id_soup = BeautifulSoup(rData, 'xml')

    if id_soup is None:
        print("No taxonomy xml returned")
        return

    if id_soup.eSearchResult is None or id_soup.eSearchResult.IdList.Id is None:
        print("No id in taxonomy xml for: " + str(rank_value) + " " + str(previous_rank_value))
        return

    id = id_soup.eSearchResult.IdList.Id.string

    #####################################################
    # Get taxonomy xml for the id
    #####################################################
    # http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&id=34886&retmode=xml
    url = '/entrez/eutils/efetch.fcgi?db=taxonomy&id=' + id + '&retmode=xml'
    response = get_response(url)
    rData = response.read()
    fastaSoup = BeautifulSoup(rData, 'xml')

    ###################################################
    # check that the genus and species names are correct
    ###################################################

    orgList = fastaSoup.find("ScientificName")  # start value in list format

    organism = orgList.contents[0]  # Genus species

    org_name = str(organism).split()

    Taxon = fastaSoup.find_all("Taxon")

    if org_name[0] != rank_value and org_name[0] != previous_rank_value:
        correct_organism = False
        Synonyms = Taxon[0].find_all("Synonym")
        Synonyms += Taxon[0].find_all("GenbankSynonym")

        for synonym in Synonyms:
            if synonym is None or synonym.string is None:
                continue

            scientific_name = rank_value + " " + previous_rank_value
            print( "scientific_name is " + str(scientific_name) + " " + ", synonym is " + str(synonym))
            if scientific_name in synonym.string:
                correct_organism = True
                break
            #elif rank_value in synonym.string:
            #        correct_organism = True
            #       break


        if correct_organism == False:
            print("synonym not found for " + str(rank_value) + " " + str(previous_rank_value))
            return

    ###############################################
    # Get the lineage
    ###############################################

    for taxa in Taxon:

        rank = taxa.find("Rank")
        # if rank.string == "no rank":
        # AllOtherRank += [str(taxa.ScientificName.string)]
        if rank.string == "superkingdom":
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
        elif rank.string == "tribe":
            Tribe = taxa.ScientificName.string
        elif rank.string == "genus":
            continue
        elif rank.string == "species":
            continue
        # Subspecies is taken from the organism name in the xml, usually None is given
        # else:
        elif rank.string != "no rank":
            print("unknown rank " + str(rank.string) + " " + str(taxa.ScientificName.string))
            # AllOtherRank += [str(taxa.ScientificName.string)]

            # if rank.string=="no rank" and taxa.ScientificName.string=="Cetartiodactyla":
            #    order = "Cetartiodactyla"

    return OrderedDict([
        ("Tribe", Tribe),
        ("Subfamily", Subfamily),
        ("Family", Family),
        ("Superfamily", Superfamily),
        ("Parvorder", Parvorder),
        ("Infraorder", Infraorder),
        ("Suborder", Suborder),
        # ("AllOtherRank", str(AllOtherRank)),
        ("Order", Order),
        ("Superorder", Superorder),
        ("Sublass", Subclass),
        ("Class", Class),
        ("Subphylum", Subphylum),
        ("Phylum", Phylum),
        ("SuperPhylum", Superphylum),
        ("Kingdom", Kingdom),
        ("Superkingdom", Superkingdom)])


def checkDb(rank, rankValue, geneList):
    print("checking database for " + rankValue)
    missing_gene_list = []
    cur = dbConnection.cursor()
    fasta = ""

    for gene in geneList:
        # check if the gene is present at the given taxonomic level
        if rank == "Species":  # look for genus and species, since two Genus can potentially have same species
            select = 'SELECT * FROM RefGenes WHERE Genus=? AND Species=? AND GeneName=?'
            cur.execute(select, (genus, species, gene))
        else:
            # select = 'SELECT * FROM RefGenes WHERE Genus=? AND GeneName=?'
            select = 'SELECT * FROM RefGenes WHERE "' + rank + '"=? AND GeneName=?'
            cur.execute(select, (rankValue, gene))
        results = cur.fetchone()

        if results is not None:  # returns one row, assumes cellCount=24 but this will change if attributes change
            # print("fetching results for " + str(results["Genus"]) + " " + str(results["Species"]) + "for gene " + gene)

            sequence = str(results["GeneSequence"])
            acc_num = str(results["AccessionNumber"])
            genus_used = str(results["Genus"])
            species_used = str(results["Species"])
            partial_flag = str(results["PartialFlag"])
            if partial_flag == "0":
                partial_or_complete = "complete"
            if partial_flag == "1":
                partial_or_complete = "partial"
            # print("fasta-ing")
            header = ">" + gene + '_' + genus_used + '_' + species_used + '_' + acc_num + '_' + rank + '_' + \
                     rankValue + '_' + partial_or_complete
            fasta += header + "\n" + sequence + "\n"

        else:  # e.g., no gene present for organism of the given rank
            missing_gene_list += [gene]  # still missing, keep in list
            # print("making the missing")
    # print("geneList is " + str(geneList) + "and missing_gene_list is " + str(missing_gene_list))
    return missing_gene_list, fasta


if __name__ == '__main__':
    # Add ability to limit the rank that a relative seed can be extracted.

    single_seed = "n"
    i = 0
    try:
        sample_list = sys.argv[1]  # raw data list from Alan's lab [id, family, genus, species]


    except IndexError:
        # print("please specify a filename")

        # sample_list = "Zaher18Samples.csv"
        sample_list = "SampleList.csv"
        print("using default sampleList " + str(sample_list))

    try:
        single_seed = sys.argv[2]
        if single_seed == "y" or single_seed == "Y":
            single_seed = "y"

    except IndexError:
        print("Gene seeds will be contained in the same file.")

    try:
        geneList = literal_eval(sys.argv[3])  # i.e. geneList=[\'COX1\', \'ND2\',\'12S\',\'16S\',\'ND5\']
        # populateRelatives = sys.argv[3]
        # Check that genes are qualified geneList names, to avoid redundant DB entries
        for gene in geneList:
            if gene.upper() == "COX1":
                continue
            if gene.upper() == "ND2":
                continue
            if gene.upper() == "12S":
                continue
            if gene.upper() == "16S":
                continue
            if gene.upper() == "ND5":
                continue
            if gene.upper() == "CYTB":
                continue
            if gene.upper()=="matK":
                continue
            if gene.upper()=="rbcL":
                continue
            # if gene.upper()=="";
            #    continue
            else:
                print("Note that only following gene names have synonyms available: ['COX1', 'ND2','12S','16S','ND5', matK, rbcL]")

    except IndexError:
        print("Using default gene list")
        #geneList = ['rbcL','matK', 'trnH', 'atpF', 'psbA', 'psbK', 'rpl32', 'rpoC1', 'rpoB', 'atpB', 'ndhF', 'rps16', 'trnC', 'trnE', 'trnG', 'trnK', 'trnS','trnT', 'trnY', 'ycf1', 'ycf6', 'COX1', 'ITS1','ITS2']
        #geneList = ['COX1', 'CYTB', 'ND2', '12S', '16S', 'ND5']
        #Plant Gene List
        geneList = ['atpB', 'atpB+rbcL', 'atpE', 'matK', 'ndhF', 'psbA+trnH',
                    'psbB', 'psbH', 'psbN', 'psbZ', 'psbZ+trnG',
                    'rbcL', 'rpl16', 'rpl32+trnL', 'rpoC1', 'rps4', 'rps16', 'trnD+trnY', 'trnK+matK',
                    'trnL+trnF', 'trnT+trnL', 'trnS+psbZ', 'ycf1', 'ITS1', 'ITS2']
        '''
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
        # populateRelatives = "No"
        '''

    try:
        max_rank = sys.argv[4]

    except IndexError:
        max_rank = 'Order'
        print("max rank is set to default: Order")

    with open(sample_list, 'r') as sample:
        lines = sample.readlines()

        for line in lines:
            # read in list of Genus species names
            # Get the genus and species name from comma delimited file (e.x. ZaherRaw.csv)
            # line = line.strip() #should not have any spaces
            fasta_lines = None
            line = line.replace('\n', '')
            cells = line.split(',')
            id = cells[0]
            family = cells[1]
            genus = cells[2]
            species = cells[3]
            trans_genus = None
            trans_species = None
            fasta = ""
            missing_gene_list = list(geneList)  # LISTS ARE MUTABLE



            ##### To do: Look for same subspecies in sqlite first ###

            if species != "sp.":
                missing_gene_list, fasta_lines = checkDb("Species", species, missing_gene_list)

                if len(fasta_lines) != 0:
                    print("checkDB for species returned " + fasta_lines)
                    print("for species " + species)
                    fasta += fasta_lines  # Two lines appended, a header and a gene sequence (may be None)
                else:
                    # Check synonyms in DB, if exists then use that name
                    trans_genus, trans_species = checkSynonyms(genus, species)
                    if trans_species is not None and trans_species!=species:
                        print(" translated " + str(genus) + " to " + str(trans_genus))
                        genus = trans_genus
                        print(" translated " + str(species) + " to" + str(trans_species))
                        species = trans_species
                        missing_gene_list, fasta_lines = checkDb("Species", species, missing_gene_list)

                        if len(fasta_lines) != 0:
                            print("checkDB for trans_genus returned " + fasta_lines)
                            print("for species " + species)
                            fasta += fasta_lines
                        else:
                            print("################# no value found in database for species " + species)

            if len(missing_gene_list) != 0:
                print("checking genus")
                missing_gene_list, fasta_lines = checkDb("Genus", genus, missing_gene_list)
                if len(fasta_lines) != 0:
                    print("checkDB returned " + fasta_lines)
                    print("for genus " + genus)
                    fasta += fasta_lines  # Two lines appended, a header and a gene sequence (may be None)
                else:
                    print("################# no value found in database for genus " + genus)

            if len(missing_gene_list) != 0:  # No genus present, search for full taxonomy
                # If genes remain to be found, find taxonomy and use relative for those genes
                # AllOtherRank,Superkingdom,Kingdom,Superphylum,Phylum,Subphylum,Class,Superorder,Order,Suborder,Infraorder,Parvorder,Superfamily,Family,Subfamily = get_full_taxonomy(genus,species)
                # ranks = [Subfamily, Family, Superfamily, Parvorder, Infraorder, Order, AllOtherRank, Superorder, Class, Subphylum, Phylum, SuperPhylum, Kingdom, SuperKingdom]
                temp_taxonomy = None
                taxonomy = OrderedDict([("Genus", genus)])  # add genus

                temp_taxonomy = get_taxonomy(genus, species)

                if temp_taxonomy is None:
                    print("####### No taxonomy for " + str(genus) + " " + str(species))
                    temp_taxonomy = get_taxonomy(family, genus)

                    if temp_taxonomy is None:
                        print("###########No taxonomy for " + str(family) + " " + str(genus))
                        temp_taxonomy = get_taxonomy(family, "")

                        if temp_taxonomy is None:
                            print("################No taxonomy for " + str(family) )
                            print("Can not search above family level for " + family + " " + genus + " " + species)
                            temp_taxonomy = OrderedDict([("Family", family)])



                taxonomy.update(temp_taxonomy)  # add everything else

                print("Taxonomy is " + str(taxonomy))

                for rank, rankValue in taxonomy.items():

                    if rank == max_rank:
                        print "max rank reached, no result for gene " + missing_gene_list[0]
                        print "missing gene list was " + missing_gene_list
                        missing_gene_list = missing_gene_list[1:]
                        print "missing gene list is now " + missing_gene_list


                    if rankValue is None:
                        continue
                    else:
                        missing_gene_list, fasta_lines = checkDb(rank, rankValue, missing_gene_list)
                        if len(fasta_lines) != 0:
                            print("length of fasta_lines is " + str(len(fasta_lines)))
                            print("checkDB returned " + fasta_lines)
                            print("for rank " + rank + " and rankValue " + rankValue)
                            fasta += fasta_lines  # Two lines appended, a header and a gene sequence (may be None)

                        if len(missing_gene_list) == 0:
                            # save (append) gene in fasta file for that species
                            # line of info followed by line of sequence
                            # filePrefix = genus + '_' + species #Want to use id instead
                            # f_out= open(filePrefix + ".seeds", 'w')
                            if single_seed == "y":
                                f_out = open("./seeds/" + id + "_" + str(i) + ".seeds", 'w')
                                i = i + 1
                            else:
                                f_out = open("./seeds/" + id + ".seeds", 'w')
                            f_out.write(fasta)
                            f_out.close()
                            break


            else:  # missing_gene_list is None, species or Genus was present
                # print("writing to file")
                # filePrefix = genus + '_' + species
                print "No genes Left"
                if single_seed == "y":
                    f_out = open("./seeds/" + id + "_" + str(i) + ".seeds", 'w')
                    i = i + 1
                else:
                    f_out = open("./seeds/" + id + ".seeds", 'w')
                f_out.write(fasta)
                f_out.close()
                #### NOTE: Still need to add an optional check to NCBI when not found in database ####
                ##### Will write over fasta files of same species with a different subspecies
                #### Need to handle AllOtherRank for optional NCBI search

        dbConnection.close()


        # Create an option to search for relatives in NCBI and populate database if none found
        # modify superextender to use only the line of the gene you want (sys args)
