import json, sys
import kyotocabinet as kc

def build_node_hash(dmp_file, names_kch, taxdb= 'nodes.kch'):
    """
    Creates a hashed interpretation of the nodes.dmp file from NCBI using kyotocabinet.
    Requires a hashed names.dmp file be generated using build_name_hash() first
    """
    namedb= kc.DB()
    if not namedb.open(names_kch, kc.DB.OREADER):
        raise IOError
    db= TaxDB()
    if not db.open(taxdb, TaxDB.OWRITER | TaxDB.OCREATE):
        sys.stderr.write("opening error: " + str(db.error()))
    try:
        for line in dmp_file:
            line= line.strip().split('\t|\t')
            key= line[0]
            vals= {
                    'tax_id': line[0],
                    'parent_id': line[1],
                    'rank': line[2],
                    'embl_code': line[3],
                    'division_id': line[4],
                    'comments': line[-1]
                    }
            vals['name']= namedb.get(key)
            db.set(key, vals)
    finally:
        db.close()

def build_name_hash(dmp_file, namedb= 'names.kch'):
    """
    Creates a hashed interpretation of the names.dmp file from NCBI using kyotocabinet.
    Only stores scientific names, as all others are pretty much useless except for legacy reasons.
    """
    db= kc.DB()
    if not db.open(namedb, db.OWRITER | db.OCREATE):
        sys.stderr.write("opening error: " + str(db.error()))
    try:
        for line in dmp_file:
            if "scientific name" not in line:
                continue
            line= line.strip().split('\t|\t')
            key= line[0]
            name= line[1]
            db.set(key, name)
    finally:
        db.close()

def build_reverse_name_hash(dmp_file, namedb= 'reverse_names.kch'):
    """
    Creates a reverse-hash interpretation of the names.dmp file from NCBI using kyotocabinet.
    Only stores scientific names, as all others are pretty much useless except for legacy reasons.
    """
    db= kc.DB()
    if not db.open(namedb, db.OWRITER | db.OCREATE):
        sys.stderr.write("opening error: " + str(db.error()))
    try:
        for line in dmp_file:
            if "scientific name" not in line:
                continue
            line= line.strip().split('\t|\t')
            key= line[0]
            name= line[1]
            db.set(name, key)
    finally:
        db.close()

def taxonomy(tax_id, db):
    """
    Returns the taxonomy "string" of the tax_id as a list structure back to (but not including) root
    """
    node= db.get(tax_id)
    if not node:
        return ''
    tax_list= []
    while node['parent_id'] != '1':
        tax_list.insert(0, node['name'])
        node= db.get(node['parent_id'])
        if not node:
            break
    return tax_list

def search(tax_string, reverse_db):
    """
    Returns taxonomy id based on a search string
    """
    node= db.get(tax_string)
    if not node:
        return ''
    return node

def open_taxdb(kch_name, method='r'):
    """
    Returns a TaxDB with safe opening
    """
    db= TaxDB()    
    if method == 'w':
        if not db.open(kch_name, db.OWRITER | db.OCREATE):
            sys.stderr.write("opening error: " + str(db.error()))
        return db
    elif method == 'r':
        if not db.open(kch_name, db.OREADER):
            sys.stderr.write("opening error: " + str(db.error())) 
        return db
    else:
        raise ValueError

class TaxDB(kc.DB):
    """
    TaxDB extends the KyotoCabinet database for storage and retrieval of complex data-structures using json
    """

    def get(self, key):
        """
        Retrieve value by key with interpretation of value as json
        """
        result= super(self.__class__, self).get(key)
        if result:
            result= json.loads(result)
        return result

    def set(self, key, value):
        """
        Store value by key with interpretation of value as json
        """
        value= json.dumps(value)
        super(self.__class__, self).set(key, value)


def check_arguments(args):
    if args.build_names and (args.tax_string or args.tax_list or args.build_nodes):
        raise argparse.error('Please create the requisite database files before specifying taxonomy IDs')
    if args.build_nodes and (args.tax_string or args.tax_list or args.build_names):
        raise argparse.error('Please create the requisite database files before specifying taxonomy IDs')
    if args.tax_string and args.tax_list:
        raise argparse.error('Specify either a single string or filename, not both')

if __name__ == '__main__':
    import argparse, os
    parser = argparse.ArgumentParser(description='Utility functions for taxonomy storage and retrieval from NCBI taxdump files', formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

        #Optional arguments
    optional = parser.add_argument_group('OPTIONAL')
    optional.add_argument('-h', action="help", help="show this help message and exit")
    optional.add_argument('-t', '--tax_string', help= 'Retreive a single taxonomy string and print to STDOUT', type=str)
    optional.add_argument('-l', '--tax_list', help= 'Process a batch of taxonomy strings read from file', type=str)
    optional.add_argument('--names', help= 'A names.kch file built using this module', type=str)
    optional.add_argument('--nodes', help= 'A nodes.kch file built using this module', type=str, required=True)
    optional.add_argument('--build_names', help= 'Build a kch for the names.dmp file', type=str)
    optional.add_argument('--build_reverse_names', help= 'Build a reverse kch for the names.dmp file', type=str)
    optional.add_argument('--build_nodes', help= 'Build a kch for the nodes.dmp file (--build_names must be run first!)', type=str)

    args = parser.parse_args()


    check_arguments(args)

    if args.build_names:
        base= os.path.splitext(args.build_names)[0]
        build_name_hash(open(args.build_names), base+'.kch')
    elif args.build_nodes:
        base= os.path.splitext(args.build_nodes)[0]
        build_node_hash(open(args.build_nodes), base+'.kch')
    elif args.build_reverse_names:
        base= os.path.splitext(args.build_reverse_names)[0]
        build_reverse_name_hash(open(args.build_reverse_names))#, base+'.kch')
    
    if args.tax_string or args.tax_list:
        taxdb= open_taxdb(args.nodes)
        if args.tax_string:
            print ';'.join(taxonomy(args.tax_string, taxdb))
        else:
            for line in open(args.tax_list):
                line= line.strip()
                print ';'.join(taxonomy(line, taxdb))
