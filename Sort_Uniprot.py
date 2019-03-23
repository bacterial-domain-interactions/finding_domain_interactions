import csv
# csvWriter method which takes in a dictionary and the name of the file you want to write it to.
# It will go through each key in your map and see if it the value is empty or not. If it is not empty then we will add
# it to the csv file you named, otherwise it will skip over to the next key.


def csvWriter(domainMap, fileName):
    """
    ---- Description
    :param domainMap:
    :param fileName:
    :return:
    """
    with open(fileName, 'w', newline='') as csvfile:
        fieldnames = ['protein', 'domain']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for key in domainMap.keys():
            if domainMap.get(key) == "Empty":
                continue
            else:
                writer.writerow({'protein': key, 'domain': domainMap.get(key)})


def give_me_my_domain_list(peter_list, which_protein, uniprot_list, domain_ID):
    '''
    ----Description of function ------
    :param peter_list:
    :param which_protein: protein 1/protein 2
    :param uniprot_list: pfam list
    :param domain_ID: pFAM/interpro
    :return:
    '''
    pFam_map = dict()
    interpro_map = dict()
    peter_set = set()
    with open(peter_list) as scanner:
        next(scanner)
        for line in scanner:
            line_split = line.split(",")
            if "1" in which_protein:
                peter_set.add(line_split[0].strip())
            elif "2" in which_protein:
                peter_set.add(line_split[1].strip())

    with open(uniprot_list) as scan:
        next(scan)
        for line in scan:
            line = line.strip()
            line = line.replace(";", " ")
            fams = line.split(",")
            i = fams[0].strip()
            for protein in peter_set:
                if i == protein.strip() or protein.strip() in i:
                    if len(fams) == 1:
                        pFam_map[fams[0]] = "Empty"
                        interpro_map[fams[0]] = "Empty"
                    if len(fams) == 2:
                        if fams[1].startswith("P"):
                            pFam_map[i] = fams[1]
                            interpro_map[i] = "Empty"
                        elif fams[1].startswith("I"):
                            pFam_map[i] = "Empty"
                            interpro_map[i] = fams[1]
                    if len(fams) == 3:
                        if fams[1].startswith("P"):
                            pFam_map[i]=fams[1]
                        if fams[1].startswith("I"):
                            interpro_map[i]=fams[2]
                        if fams[2].startswith("P"):
                            pFam_map[i] = fams[1]
                        if fams[2].startswith("I"):
                            interpro_map[i] = fams[2]
                    if fams[1] is None or fams[1] == "":
                        pFam_map[i] = "Empty"
                    if fams[2] is None or fams[2] == "":
                        interpro_map[i] = "Empty"

    domain_ID = domain_ID.upper()
    fileName= domain_ID + "_" + which_protein + ".csv"  # Should read PFAM_1.csv.

    if domain_ID == "PFAM":
        csvWriter(pFam_map, fileName)

    elif domain_ID == "INTERPRO":
        csvWriter(interpro_map, fileName)

    return fileName