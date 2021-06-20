/*************************************************************************/
/*structures.h                                                           */
/*This code contains the structures to save important results            */
/*************************************************************************/
#ifndef STRUCTURES_H
#define STRUCTURES_H

typedef struct esl
{
  double volt;
  double time;
  struct esl *sig;
}ESL;

typedef struct times
{
  double tf; // Type B
  double tr; // Type B
  double b; // Type B
  double trstartind; // Type B
  double trendind; // Type B
  double tfstartind; // Type B
  double tfendind; // Type B

  double t1; // Type A (all bit rates)
  double t1startind; // Type A (all bit rates)
  double t1start; // Type A (all bit rates)
  double t1endind; // Type A (all bit rates)

  double t2; // Type A
  double t2startind; // Type A
  double t2start; // Type A
  double t3; // Type A
  double t3end; // Type A
  double t3endind; // Type A
  double t4; // Type A
  double t4endind; // Type A

  double t5; // Type A (higher bit rates)
  double t5startind; // Type A (higher bit rates)
  double t6; // Type A (higher bit rates)
  double t6end; // Type A (higher bit rates)
  double t6endind; // Type A (higher bit rates)
  double a; // Type A (higher bit rates)
  double tploone; //Type A (higher bit rates)
}TIMES;

typedef struct shootreader
{
  double shootind;
  double shootind_b;
  double hf_reader;
  double hr_reader;
  double above;
  double above_b;
}SHOOTREADER;

#endif

