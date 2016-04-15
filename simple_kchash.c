#include <kclangc.h>
#include <string.h>
/* when compiling, make sure to use the libkyotocabinet (e.g. /usr/lib64/libkyotocabinet.so) */

int main(int argc, char* argv[])
{
	FILE *stream;
	char *key, *value, *line = NULL;
	size_t len = 0;
	int count = 0;
	ssize_t read;

	KCDB* db;
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
