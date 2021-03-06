#!/usr/bin/python

from gene_list_utils import *

'''
Parse the drugbank XML dump and use an HGNC dictionary to get a mapping of drug generic name to approved human gene symbol.
Suggested input:
    wget http://www.drugbank.ca/system/downloads/current/drugbank.xml.zip
    unzip drugbank.xml.zip
And before you use the drugbank.xml file, use sed or a text editor to delete 'xmlns="http:///www.drugbank.ca"' from the second line.
Otherwise you need to reference the namespace every time you access any node.
This output only includes 1) approved drugs with 2) known mechanisms of action with 3) human targets with 4) value gene symbols.
Note that annotations exist for other groups as well, but we are only addressing approved drugs here. Counts in 2015 were:
$ cat drugbank.xml | grep "<group>" | sort | uniq -c
# 1757     <group>approved</group>
# 5064     <group>experimental</group>
#  186     <group>illicit</group>
# 1231     <group>investigational</group>
#   89     <group>nutraceutical</group>
#  179     <group>withdrawn</group>
'''
def parse_drugbank(drugbank_path,hgnc,debug=False):
    drug_number = 1
    drug_gene = {} # keys will be drug generic names, values will be gene symbols
    drug_gene_action = {} # keys will be drug generic names, values will be list [gene symbol, action]
    drug_classifications = {} # keys will be drug generic names, values will be list [classification, atc_code]
    drug_categories = {} # keys will be drug generic names, values will be list of categories
    tree = ET.parse(drugbank_path)
    root = tree.getroot()
    drugs = root.findall("drug")
    for drug in drugs:
        # get its common name
        name = drug.find("name")
        if name is None:
            continue
        # figure out if approved
        is_approved = False
        groups = drug.find("groups")
        if groups is not None:
            for group in groups.findall("group"):
                if group.text == "approved":
                    is_approved = True
        if not is_approved:
            continue
        generic_name = name.text.lower()
        if debug:
            sys.stderr.write('\rdrug ' + str(drug_number) + ': ' + generic_name)
            drug_number += 1
        gene_symbol = ''
        action_text = ''
        targets_element = drug.find("targets")
        if targets_element is not None:
            targets = targets_element.findall("target")
            for target in targets:
                # not all targets have a position, but for those that do,
                # only take the top-ranked target association
                if target.attrib.has_key("position"):
                    if target.attrib["position"] != '1':
                        continue
                known_action = target.find("known-action")
                if known_action.text != "yes":
                    continue
                else:
                    actions = target.find("actions")
                    if actions is not None:
                        action = actions.find("action")
                        if action is not None:
                            action_text = action.text
                polypeptide = target.find("polypeptide")
                if polypeptide is None:
                    continue
                organism = polypeptide.find("organism")
                if organism.text != "Human":
                    continue
                extids = polypeptide.find("external-identifiers")
                for extid in extids.findall("external-identifier"):
                    if extid.find("resource").text == "HUGO Gene Nomenclature Committee (HGNC)":
                        hgnc_id = extid.find("identifier").text
                if hgnc_id is not None:
                    if hgnc.has_key(hgnc_id):
                        gene_symbol = hgnc[hgnc_id]['symbol']
                        break
            drug_gene[generic_name] = gene_symbol # add key-value pair to dictionary
            drug_gene_action[generic_name] = [gene_symbol, action_text] # add key-value pair to dictionary
            # figure out categories
            drug_categories[generic_name] = []
            cats_header = drug.findall("categories")
            if cats_header is not None:
                cats = cats_header[0].findall("category")
                for cat in cats:
                    cat_text = cat.find("category").text
                    drug_categories[generic_name].append(cat_text)
            # figure out type
            dtype = drug.attrib['type']
            # extract atc code
            atc_code_text = ''
            atc_codes_parent = drug.findall("atc-codes")
            if atc_codes_parent is not None:
                atc_codes = atc_codes_parent[0].findall("atc-code")
                if atc_codes is not None and len(atc_codes) > 0:
                    atc_code = atc_codes[0]
                    atc_code_text = atc_code.attrib['code']
            drug_classifications[generic_name] = [dtype, atc_code_text]
            if debug:
                if gene_symbol is None:
                    print ""
                else:
                    print gene_symbol
    return [drug_gene, drug_gene_action, drug_classifications, drug_categories]

if __name__ == '__main__':
    hgnc = parse_hgnc("../drug_target_lof_full/raw_data/gene_with_protein_product_2018_09_13.txt",mode='id')
    parsed = parse_drugbank("../drug_target_lof_full/raw_data/drugbank.xml",hgnc,debug=True)
    drug_gene = parsed[0]
    drug_gene_action = parsed[1]
    drug_classifications = parsed[2]
    drug_categories = parsed[3]
    with open('data/drug_gene_match.tsv',mode='w') as f:
        f.write("drug\tgene\n")
        print_dict(drug_gene,f)
    with open('data/drug_gene_action_match.tsv',mode='w') as f:
        f.write("drug\tgene\taction\n")
        print_dict_with_list(drug_gene_action,f)
    with open('data/drug_classifications.tsv',mode='w') as f:
        f.write("drug\ttype\tatc\n")
        print_dict_with_list(drug_classifications,f)
    with open('data/drug_categories.tsv',mode='w') as f:
        f.write("drug\tcategory\n")
        print_dict_relational_list(drug_categories,f)
    with open('lists/fda_approved_drug_targets.tsv',mode='w') as f:
        print_values(drug_gene,f)
