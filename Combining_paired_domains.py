import csv


def getMap(fileName):
    protein_map = dict()
    with open(fileName) as scanner:
        next(scanner)
        for line in scanner:
            domains = line.split(",")
            protein_map[domains[0]] = domains[1]
    return protein_map

# protein1_file = "PFAM_protein 1.csv"
# protein1 = getMap(protein1_file)
#
# protein2_file = "PFAM_protein 2.csv"
# protein2=getMap(protein2_file)

# paired_protein_list = "HelicoBacter_paired.csv"

def get_me_my_paired_domain_list(paired_protein_list, protein1_file, protein2_file):
    fileName = "combined.csv"
    protein1 = getMap(protein1_file)
    protein2 = getMap(protein2_file)
    with open(fileName, 'w', newline="") as csvfile, open(paired_protein_list) as scan:
        fieldname = ["ProteinX", "XDomains", "ProteinY", "YDomains"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldname)
        writer.writeheader()
        next(scan)
        for line in scan:
            proteins=line.split(",")
            if proteins[0].strip() in protein1.keys() and proteins[1].strip() in protein2.keys():
                writer.writerow({'ProteinX': proteins[0].strip(),'XDomains': protein1.get(proteins[0].strip()).strip(),
                                 'ProteinY': proteins[1].strip(), 'YDomains': protein2.get(proteins[1].strip()).strip()})
    return fileName