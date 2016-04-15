/*import kyotocabinet as kc
import sys
import argparse, os

parser = argparse.ArgumentParser(description='Simple generation of kyotocabinet hash table from tab-separated table', 
        formatter_class=argparse.ArgumentDefaultsHelpFormatter, add_help=False)

#Required arguments
required = parser.add_argument_group('REQUIRED')
required.add_argument('-t', help= 'input tsv', required=True, type=str)
required.add_argument('-k', help= 'output kyotocabinet', required=True, type=str)

#Optional arguments
optional = parser.add_argument_group('OPTIONAL')
optional.add_argument('-h', action="help", help="show this help message and exit")

args = parser.parse_args()

db= kc.DB()
if not db.open(args.k, db.OWRITER | db.OCREATE):
    sys.stderr.write("opening error: " + str(db.error()))

for line in open(args.t):
    key, value= line.split('\t')[:2]
    value= value.strip()
    db.set(key, value)

if not db.close():
    sys.stderr.write("close error: " + str(db.error()))
    */
#include <kclangc.h>
#include <string.h>

int main(int argc, char* argv[])
{
	FILE *stream;
	char *key, *value, *line = NULL;
	size_t len = 0;
	int count = 0;
	ssize_t read;
	KCDB* db;
	/*
	if (stream == NULL)
		exit(EXIT_FAILURE);

	while ((read = getline(&line, &len, stream)) != -1) {
		key = strsep(&line, "\t");
		value = strsep(&line, "\t");
		printf("Retrieved line of length %zu :\n", read);
		printf("%s\t", key);
		printf("%s\n", value);
	}

	int32_t t = kcdbopen(db, argv[2], KCOWRITER|KCOCREATE);
	kcdbclose(db);
	*/
	db = kcdbnew();
	if (!kcdbopen(db, argv[2], KCOWRITER | KCOCREATE))
		fprintf(stderr, "open error: %s\n", kcecodename(kcdbecode(db)));

	stream = fopen(argv[1], "r");
	if (stream == NULL)
		exit(EXIT_FAILURE);

	while ((read = getline(&line, &len, stream)) != -1) {
		key = strsep(&line, "\t");
		value = strsep(&line, "\t");
		if (!kcdbset(db, key, strlen(key), value, strlen(value)))
			fprintf(stderr, "set error: %s\n", kcecodename(kcdbecode(db)));
		count++;
		if (count%1000000 == 0)
			fprintf(stderr, "%d lines processed\n", count);
	}

	free(line);
	fclose(stream);
	/* close the database */
	if (!kcdbclose(db)) {
		fprintf(stderr, "close error: %s\n", kcecodename(kcdbecode(db)));
	}
	kcdbdel(db);

}
