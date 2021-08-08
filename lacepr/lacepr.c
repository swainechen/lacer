#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <float.h>
#include <math.h>
#include <ctype.h>
#include "sam.h"
#include "sam_header.h"
#include "htslib/sam.h"
#include "lacepr.h"

const char *version = "0.11";

cache_t *cache;

void init_cache (int n) {
  int i, j, k, l;
  cache = realloc(cache, n * sizeof(cache_t));
  for (i=0; i < n; i++) {
    for (j=0; j < MAX_REASONABLE_Q+1; j++) {
      for (k=0; k < MAX_CONTEXT+1; k++) {
        for (l=0; l < MAX_CYCLE+1; l++) {
          cache[i][j][k][l] = 0;
        }
      }
    }
  }
}

long choose (int n, int k) {
  long result = 1;
  int j = 1;
  if (k == n || k == 0) {
    return 1;
  }
  if (k > n || k < 0 || n < 1) {
    return 0;
  }
  if (n - k < k) {
    k = n - k;
  }
  while (j < k) {
    result *= n/j;
    n--;
    j++;
  }
}

double log10binomial (double kd, int n, double p) {
  int k;
  double mean;
  double stdev;
  k = lround(kd);
//printf("kd %f, n %d, p %f\n", kd, n, p);
  if (k >= 0 && n > 0 && p > 0 && p < 1) {
    if (n > 100) {
      // use normal approximation
      mean = n * p;
      stdev = sqrt(mean * (1 - p));
      return (-0.5 * ( (k - mean) * (k - mean) / n / p / (1-p) ) / log(10) - 0.5 * log10(2 * M_PI * n * p * (1-p)));
    } else {
      return ( log10(choose(n, k)) + k * log10(p) + (n-k) * log10(1-p) );
    }
  }
  return 0;
}

int delta (int origq, int obs, double err) {
  int i;
  int maxi = -1;
  double temp;
  double max = -DBL_MAX;
  double p[MAX_REASONABLE_Q + 1];
  for (i = MIN_Q; i <= MAX_REASONABLE_Q; i++) {
    if (abs(origq - i) > MAX_USABLE_Q) {
      p[i] = prior[MAX_USABLE_Q];
    } else {
      p[i] = prior[abs(origq - i)];
    }
    temp = log10binomial(err, obs, pow(10, -i/10.0));
    if (temp != 0) {
      p[i] += temp;
    }
    if (p[i] > max) {
      maxi = i;
      max = p[i];
    }
//printf ("origq %d, i %d, diff %f temp %f pi %f\n", origq, i, prior[abs(origq - i)], temp, p[i]);
  }
  if (maxi >= MIN_Q && maxi <= MAX_REASONABLE_Q) {
    return (maxi - origq);
  } else {
    return 0;
  }
}

int newq (int8_t origq, char *rg, int cycle, char preceding, char current, recal_t *recaldata, char** rglist, int num_rg) {
  int nq = 40;
  int running_q;
  int adjust_q = 0;
  int adjust_cyc = 0;
  int adjust_con = 0;
  char context[3];
  context[0] = preceding;
  context[1] = current;
  context[2] = '\0';
  int con_i = get_context_index(context);
  int cyc_i = get_cycle_index(cycle);

  int rgindex = get_rg_index(rglist, num_rg, rg);
  if (rgindex == -1) {
    return origq;
  }

  if (cache[rgindex][origq][con_i][cyc_i] > 0) {
    return cache[rgindex][origq][con_i][cyc_i];
  }

  running_q = recaldata[rgindex].Quality;
  if (recaldata[rgindex].OrigQual[origq].Observations > 0) {
    adjust_q = delta(running_q, recaldata[rgindex].OrigQual[origq].Observations, recaldata[rgindex].OrigQual[origq].Errors);
//printf("base %d, obs %d, err %f, adjust_q %d\n", running_q, recaldata[rgindex].OrigQual[origq].Observations, recaldata[rgindex].OrigQual[origq].Errors, adjust_q);
  }
  running_q += adjust_q;
 
  if (recaldata[rgindex].Context[con_i].OrigQual[origq].Observations > 0) {
    adjust_con = delta(running_q, recaldata[rgindex].Context[con_i].OrigQual[origq].Observations, recaldata[rgindex].Context[con_i].OrigQual[origq].Errors);
//printf("base %d, context %s, index %d, obs %d, err %f, adjust_q %d\n", running_q, context, con_i, recaldata[rgindex].Context[con_i].OrigQual[origq].Observations, recaldata[rgindex].Context[con_i].OrigQual[origq].Errors, adjust_con);
  }
  if (recaldata[rgindex].Cycle[get_cycle_index(cycle)].OrigQual[origq].Observations > 0) {
    adjust_cyc = delta(running_q, recaldata[rgindex].Cycle[cyc_i].OrigQual[origq].Observations, recaldata[rgindex].Cycle[cyc_i].OrigQual[origq].Errors);
//printf("base %d, cycle %d, index %d, obs %d, err %f, adjust_q %d\n", running_q, cycle, cyc_i, recaldata[rgindex].Cycle[cyc_i].OrigQual[origq].Observations, recaldata[rgindex].Cycle[cyc_i].OrigQual[origq].Errors, adjust_cyc);
  }
  nq = running_q + adjust_con + adjust_cyc;
//  printf ("orig %d prior %d, cycle %d, context %s, adjust %d %d %d, new %d\n\n\n", origq, recaldata[rgindex].Quality, cycle, context, adjust_q, adjust_con, adjust_cyc, nq);
  if (nq == 0) {
    nq = origq;
  }
  if (nq < MIN_Q) {
    nq = MIN_Q;
  }
  if (nq > MAX_REASONABLE_Q) {
    nq = MAX_REASONABLE_Q;
  }
  cache[rgindex][origq][con_i][cyc_i] = nq;
  return nq;
}

int get_rg_index (char** rglist, int num_rg, char* rg) {
  int i;
  for (i=0; i < num_rg; i++) {
    if (strcmp(rg, rglist[i]) == 0) {
      return(i);
    }
  }
  return(-1);
}

// we're really only going to assume there is 2 base context
int get_context_index (char *s) {
  int i = 0;
  int sum = 0;
  while (s[i] != '\0') {
    sum *= 4;
    if (s[i] == 'A') {
      // this is 0
    } else if (s[i] == 'C') {
      sum += 1;
    } else if (s[i] == 'G') {
      sum += 2;
    } else if (s[i] == 'T') {
      sum += 3;
    } else {
      return 0;	// bad data
    }
    i++;
  }
  sum++;	// we want to give 1-16, not 0-15
  return sum;
}

// cycle here should be an integer, positive or negative
// encode negatives as even numbers
// 0 is bad data
int get_cycle_index (int i) {
  if (i > 0 && i < MAX_CYCLE) {
    return (i*2-1);
  } else if (i < 0 && -i < MAX_CYCLE) {
    return (-i*2);
  }
  return 0;	// bad data
}

// assume we just get an unallocated pointer, make memory and initialize
recal_t* init_recal (int num_rg) {
  int i, j, k;
  recal_t *data = malloc(sizeof(recal_t) * num_rg);
  for (i=0; i < num_rg; i++) {
    data[i].Quality = 0;
    data[i].Observations = 0;
    data[i].Errors = 0;
    for (j=0; j < MAX_Q+1; j++) {
      data[i].OrigQual[j].Quality = 0;
      data[i].OrigQual[j].Observations = 0;
      data[i].OrigQual[j].Errors = 0;
    }
    for (j=0; j < MAX_CONTEXT; j++) {
      for (k=0; k < MAX_Q+1; k++) {
        data[i].Context[j].OrigQual[k].Quality = 0;
        data[i].Context[j].OrigQual[k].Observations = 0;
        data[i].Context[j].OrigQual[k].Errors = 0;
      }
    }
    for (j=0; j < MAX_CYCLE; j++) {
      for (k=0; k < MAX_Q+1; k++) {
        data[i].Cycle[j].OrigQual[k].Quality = 0;
        data[i].Cycle[j].OrigQual[k].Observations = 0;
        data[i].Cycle[j].OrigQual[k].Errors = 0;
      }
    }
  }
  return data;
}

void read_recal (char* file, char** rglist, int num_rg, recal_t *data) {
  char *in = 0;
  size_t n = 0;
  ssize_t bytes;
  char *tkn;
  int i;
  int parsed;
  char rg[MAX_FIELD];
  char covariate[MAX_FIELD];
  char covname[MAX_FIELD];
  double empqual;
  double estqual;
  int qual;
  char eventtype[4];
  int obs;
  double err;
  int rgindex;
  int covindex;
  FILE *r = fopen(file, "r");

  while ( (bytes = getline(&in, &n, r)) != -1 ) {
    if (bytes > 0) {
      tkn = strtok(in, ":");
      if (tkn[0] == '#') {
        while (tkn != NULL) {
          tkn = strtok(NULL, ":");
          if (tkn != NULL) {
            if (strcmp(tkn, "RecalTable0") == 0) {
              // RecalTable0 should be:
              // 0 - ReadGroup
              // 1 - EventType
              // 2 - EmpiricalQuality
              // 3 - EstimatedQReported
              // 4 - Observations
              // 5 - Errors
              bytes = getline(&in, &n, r);
              i = 0;
              while ((bytes = getline(&in, &n, r)) != -1 && !whitespaceline(in)) {
                parsed = sscanf(in, "%s %s %lf %lf %d %lf", rg, eventtype, &empqual, &estqual, &obs, &err);
                if (parsed == 6) {
                  rgindex = get_rg_index(rglist, num_rg, rg);
                  data[rgindex].Quality = lround(empqual);
                  data[rgindex].Observations = obs;
                  data[rgindex].Errors = err;
//printf("rg %s index %d obs %d err %f\n", rg, rgindex, data[rgindex].Observations, data[rgindex].Errors);
//printf("rg %s index %d obs %d err %f\n", rg, rgindex, obs, err);
                } else {
                  printf("Parse error for RecalTable0 on line:\n%s", in);
                }
              }
            } else if (strcmp(tkn, "RecalTable1") == 0) {
              // RecalTable1 should be:
              // 0 - ReadGroup
              // 1 - QualityScore
              // 2 - EventType
              // 3 - EmpiricalQuality
              // 4 - Observations
              // 5 - Errors
              bytes = getline(&in, &n, r);
              while ((bytes = getline(&in, &n, r)) != -1 && !whitespaceline(in)) {
                parsed = sscanf(in, "%s %d %s %lf %d %lf", rg, &qual, eventtype, &empqual, &obs, &err);
                if (parsed == 6 && strcmp(eventtype,"M") == 0) {
                  rgindex = get_rg_index(rglist, num_rg, rg);
                  data[rgindex].OrigQual[qual].Quality = lround(empqual);
                  data[rgindex].OrigQual[qual].Observations = obs;
                  data[rgindex].OrigQual[qual].Errors = err;
                } else {
                  printf("Parse error for RecalTable1 on line:\n%s", in);
                }
//                printf("Table 1: %d, %s, %d, %s, %f, %d, %f\n", parsed, rg, qual, eventtype, empqual, data[rgindex].OrigQual[qual].Observations, data[rgindex].OrigQual[qual].Errors);
              }
            } else if (strcmp(tkn, "RecalTable2") == 0) {
              // RecalTable2 should be:
              // 0 - ReadGroup
              // 1 - QualityScore
              // 2 - CovariateValue
              // 3 - CovariateName
              // 4 - EventType
              // 5 - EmpiricalQuality
              // 6 - Observations
              // 7 - Errors
              bytes = getline(&in, &n, r);
              while ((bytes = getline(&in, &n, r)) != -1 && !whitespaceline(in)) {
                parsed = sscanf(in, "%s %d %s %s %s %lf %d %lf", rg, &qual, covariate, covname, eventtype, &empqual, &obs, &err);
                if (parsed == 8 && strcmp(eventtype,"M") == 0) {
                  rgindex = get_rg_index(rglist, num_rg, rg);
                  if (strcmp(covname, "Cycle") == 0) {
                    // covariate is a string but if it's cycle it's a number
                    if ((covindex = get_cycle_index(atoi(covariate))) != -1) {
                      data[rgindex].Cycle[covindex].OrigQual[qual].Quality = lround(empqual);
                      data[rgindex].Cycle[covindex].OrigQual[qual].Observations = obs;
                      data[rgindex].Cycle[covindex].OrigQual[qual].Errors = err;
//                  printf("Table 2: %d, %s, %d, %s, %s, %s, %f, %d, %f\n", parsed, rg, qual, covariate, covname, eventtype, empqual, data[rgindex].Cycle[covindex].OrigQual[qual].Observations, data[rgindex].Cycle[covindex].OrigQual[qual].Errors);
                    }
                  } else if (strcmp(covname, "Context") == 0) {
                    if ((covindex = get_context_index(covariate)) != -1) {
                      data[rgindex].Context[covindex].OrigQual[qual].Quality = lround(empqual);
                      data[rgindex].Context[covindex].OrigQual[qual].Observations = obs;
                      data[rgindex].Context[covindex].OrigQual[qual].Errors = err;
//                    printf("Table 2: %d, %s, %d, %s, %s, %s, %f, %d, %f\n", parsed, rg, qual, covariate, covname, eventtype, empqual, data[rgindex].Context[covindex].OrigQual[qual].Observations, data[rgindex].Context[covindex].OrigQual[qual].Errors);
                    }
                  }
                
                } else {
                  printf ("Parse error for RecalTable2 on line:\n%s", in);
                }
              }
            }
          }
        }
      }
    }
  }
}

int main (int argc, char *argv[])
{
  int getopt_c;
  char *inbam;
  char *outbam;
  char *recal_file;
  char *rg_field = "PU";
  while (1) {
    int option_index = 0;
    static struct option long_options[] = {
      {"field", required_argument, 0, 'f'},
      {"in",    required_argument, 0, 'i'},
      {"recal", required_argument, 0, 'r'},
      {"out",   required_argument, 0, 'o'},
      {0, 0, 0, 0}
    };
    getopt_c = getopt_long(argc, argv, "f:i:r:o:", long_options, &option_index);
    if (getopt_c == -1)
      break;
    switch (getopt_c) {
      case 'f': rg_field = optarg; break;
      case 'i': inbam = optarg; break;
      case 'r': recal_file = optarg; break;
      case 'o': outbam = optarg; break;
    }
  }
  if(strlen(rg_field) > 2) { rg_field[2] = '\0'; }

  if (access(inbam, F_OK) != 0 || access(recal_file, F_OK) != 0) {
    fprintf(stderr, "Lacepr version %s; headers %s\n", version, header_version);
    fprintf(stderr, "Usage: lacepr [ --field <PU|LB|SM> ] --in <in.bam> --recal <recal.tab> --out <out.bam>\n");
    return 1;
  }

  bam1_t *b;
  uint8_t *rawdata;
  int8_t *sequence;
  uint8_t *quality;
  uint8_t *aux;
  int bytes;
  uint8_t *rg;
  char *rg_const = "foo";
  char **rglist;
  int num_rg = 0;
  const char *key, *val;
  char *temp;
  void *iter;
  char preceding;
  char current;
  int cycle;
  int i;
  int debug = 0;
  int qlen;
  int max_length = 0;
  int max_aux = 0;
  int lines = 0;
  int8_t seq_comp_table[16] = { 0, 8, 4, 12, 2, 10, 9, 14, 1, 6, 5, 13, 3, 11, 7, 15 };
  recal_t *recaldata;
  char *q_pointer;
  uint8_t *q_temp;
// from bam_import.c
// char *bam_nt16_rev_tabl            = "=ACMGRSVTWYHKDBN";
// the appropriate complement would be
// char *bam_nt16_rev_tabl_complement = "=TGKCYSBAWRDMHVN";
// we're only going to accept "real" nts
  const char *nt16_table      = ".AC.G...T.......";
  const char *nt16_table_comp = ".TG.C...A.......";

  sequence = 0;
  quality = NULL;
  rg = 0;

  samfile_t *fp = samopen(inbam, "rb", 0);
  samfile_t *outfp = samopen(outbam, "wb", fp->header);
  bam_header_t *bh = fp->header;

  iter = sam_header_parse2(fp->header->text);
  sam_header_parse2(fp->header->text);
  rglist = malloc(MAX_RG * sizeof(char*));
  num_rg = 0;
  while (iter = sam_header2key_val(iter, "RG", "ID", rg_field, &key, &val)) {
    rglist[num_rg] = malloc(MAX_FIELD * sizeof(char));
//    strcpy(rglist[num_rg], key);
    strcpy(rglist[num_rg], val);
    num_rg++;
  }
  init_cache(num_rg);
  recaldata = init_recal(num_rg);
  read_recal(recal_file, rglist, num_rg, recaldata);
  b = bam_init1();
  while ((bytes = samread(fp, b)) > 0) {
    qlen = b->core.l_qseq;
    if (qlen + 1 > max_length) {
      max_length = qlen + 1;
      kroundup32(max_length);
      sequence = realloc(sequence, max_length);
//      quality = realloc(quality, max_length);
    }
    sequence[qlen] = '\0';
//    quality[qlen] = '\0';
    rawdata = bam_get_seq(b);
    uint8_t backup[b->l_data];
    for(i=0; i < b->l_data; i++) {
      backup[i] = b->data[i];
    }
    int offset = ((b)->core.n_cigar<<2) + (b)->core.l_qname + (((b)->core.l_qseq + 1)>>1);
    quality = backup + offset;
    rg = bam_aux_get(b, "RG");
    rg = rg + 1;
    if (!rg) {
      rg = rg_const;
    }
    for (i = 0; i < qlen; i++) {
      sequence[i] = bam_seqi(rawdata, i);
    }
    // this is how they access qualities in bam.h
    q_pointer = b->data + b->core.n_cigar * 4 + b->core.l_qname + ((b->core.l_qseq+1)>>1);
    for (i = 0; i < qlen; i++) {
      // get context and cycle
      if (b->core.flag & 16) {	// mapped to reverse strand
        current = nt16_table_comp[sequence[i]];
        if (i < qlen-1) {
          preceding = nt16_table_comp[sequence[i+1]];
        } else {
          preceding = '.';
        }
        cycle = qlen - i;
      } else {	// mapped to forward strand
        current = nt16_table[sequence[i]];
        if (i > 0) {
          preceding = nt16_table[sequence[i-1]];
        } else {
          preceding = '.';
        }
        cycle = i + 1;
      }
      // take care of cycle for first / second read
      if (b->core.flag & 128) {	// second of a pair
        cycle = -cycle;
      }
      debug && fprintf(stderr, "rg %s, cycle %d, sequence %d, pre %c, cur %c, orig %d, new %d\n", rg, cycle, sequence[i], preceding, current, quality[i], newq(quality[i], rg, cycle, preceding, current, recaldata, rglist, num_rg));
      q_pointer[i] = newq(quality[i], rg, cycle, preceding, current, recaldata, rglist, num_rg);
    }
    samwrite(outfp, b);
  }
  bam_destroy1(b);
  samclose(outfp);
  samclose(fp);
  return(0);
}
