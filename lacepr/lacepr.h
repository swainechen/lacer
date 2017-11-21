#define MAX_CYCLE 1024
#define MAX_FIELD 256
#define MAX_RG 256
#define MIN_Q 2
#define MAX_Q 93
#define MAX_USABLE_Q 40
#define MAX_REASONABLE_Q 60
#define CONTEXT_LENGTH 2
#define MAX_CONTEXT (1 << (CONTEXT_LENGTH * 2))

const char *header_version = "0.1";

const double prior[] = {
  -0.0457574905606751,	// 0
  -0.914346454367179,
  -3.52011334578669,
  -7.86305816481921,
  -13.9431809114647,
  -21.7604815857233,
  -31.3149601875948,
  -42.6066167170793,
  -55.6354511741769,
  -70.4014635588875,
  -86.904653871211,	// 10
  -105.145022111148,
  -125.122568278697,
  -146.83729237386,
  -170.289194396635,
  -195.478274347024,
  -222.404532225026,
  -251.06796803064,
  -281.468581763868,
  -313.606373424723,
  -DBL_MAX,	// 20
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,	// 30
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX,
  -DBL_MAX	// 40
};

typedef struct {
  int Quality;
  int Observations;
  double Errors;
} point_recal_t;

typedef struct {
  point_recal_t OrigQual[MAX_Q + 1];
} recal_set_t;

// we'll use an array of these - one for each read group
// then each line of the recal table will be in another array
// for example:
// recal_t *data;
// data[0] is for the 0th read group
// data[0].Observations and data[0].Errors
// data[0].OrigQual[30].Observations for table 1, reported quality 30
// data[0].Cycle[50].OrigQual[25].Observations for table 2, cycle 50, qual 25
// data[0].Context[5].OrigQual[25].Observations for table 2, context 5, qual 25
typedef struct {
  int Quality;
  int Observations;
  double Errors;
  point_recal_t OrigQual[MAX_Q + 1];
  recal_set_t Cycle[MAX_CYCLE + 1];
  recal_set_t Context[MAX_CONTEXT + 1];
} recal_t;

typedef int cache_t[MAX_REASONABLE_Q + 1][MAX_CONTEXT + 1][MAX_CYCLE + 1];

// http://stackoverflow.com/questions/3981510/getline-check-if-line-is-whitespace
int whitespaceline (const char *s) {
  while (s[0] != '\0') {
    if (!isspace( (unsigned char)*s )) return 0;
    s++;
  }
  return 1;
}

int get_rg_index (char** rglist, int num_rg, char* rg);
int get_context_index (char *s);
int get_cycle_index (int c);
