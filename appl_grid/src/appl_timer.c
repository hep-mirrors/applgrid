/**************************************************************************
 **
 **   File:         appl_timer.c        
 **
 **   Description:    
 **                   
 **                   
 ** 
 **   Author:       M.Sutton    
 **
 **   Created:      Fri Jun 29 18:59:11 BST 2001
 **   Modified:     
 **                   
 **                   
 **
 **************************************************************************/ 

#include <stdlib.h>
#include <pthread.h>

#include "appl_grid/appl_timer.h"


pthread_mutex_t time_lock = PTHREAD_MUTEX_INITIALIZER;  
static pthread_key_t   time_key;


typedef struct {
  struct timeval start_time;
  struct timeval stop_time;
  struct timeval since_time;  
  struct timeval diff_time;
} __appl_timer;



int appl_timersub(struct timeval* , struct timeval* , struct timeval* );

/** destructor **/
void __free_timer(__appl_timer* t) {  if ( t!=NULL ) free(t); }


/** initialiser **/
void appl_timer_initkey(void) {  pthread_key_create(&time_key, (void*)__free_timer); }


/** initialise the appl_timer thread specific key **/
__appl_timer* init_appl_timer(void) 
{
  __appl_timer* t = (__appl_timer*) malloc(1*sizeof(__appl_timer));
  pthread_setspecific(time_key, (void*)t);
  return t;
}


/**  get the time, (probably needlessly) mutex protected 
 **  call to gettimeofday
 **/

void __appl__gettime(struct timeval* t)
{
  pthread_mutex_lock(&time_lock);
  gettimeofday (t, NULL);            
  pthread_mutex_unlock(&time_lock);
}




/**  store the times in thread specific keys to avoid conflicts 
 **  of the timer routines from different threads 
 **/

double appl_timer(int flag)
{
  double      time = 0;
  __appl_timer* t = (__appl_timer*)pthread_getspecific(time_key);
  
  if (flag==APPL_TIMER_START) {
    if (t==NULL) t = init_appl_timer();
    __appl__gettime (&(t->start_time));
    t->since_time = t->start_time;
  }
  if (flag==APPL_TIMER_STOP)  {
    if (t==NULL) return -1;  /** if the timer wasn't initialised ... **/
    __appl__gettime (&(t->stop_time));
    t->since_time = t->stop_time;
    appl_timersub( &(t->stop_time), &(t->start_time), &(t->diff_time) );  
    time = (t->diff_time.tv_sec*1000.0) + (t->diff_time.tv_usec/1000.0);
  }
  if (flag==APPL_TIMER_SINCE)  {
    if (t==NULL) return -1;  /** if the timer wasn't initialised ... **/
    __appl__gettime (&(t->stop_time));
    appl_timersub( &(t->stop_time), &(t->since_time), &(t->diff_time) );  
    t->since_time = t->stop_time;
    time = (t->diff_time.tv_sec*1000.0) + (t->diff_time.tv_usec/1000.0);
  }
  
  return time;
}




/**  significantly (0.02ms) faster and simpler timer start/stop 
 **  functions 
 **/

struct timeval appl_timer_start(void)
{
  struct timeval start_time;
  __appl__gettime (&start_time);
  return start_time;
}


double appl_d(struct timeval time) { return (time.tv_sec*1000.0) + (time.tv_usec/1000.0); }

double appl_timer_stop(struct timeval start_time)
{
  double time = 0;
  struct timeval stop_time;
  struct timeval diff_time;
  double diff;

  __appl__gettime (&stop_time);
  pthread_mutex_lock(&time_lock);
  appl_timersub( &stop_time, &start_time, &diff_time );
  //  printf("diff: %d\n", _timersub( &stop_time, &start_time, &diff_time ) );
  //  diff = appl_d(stop_time)-appl_d(start_time);
  //  printf("timer:   %12.5f %12.5f  %12.5f   %12.5f\n", appl_d(stop_time), appl_d(start_time), appl_d(diff_time), diff ); 
  pthread_mutex_unlock(&time_lock);  
  time = (diff_time.tv_sec*1000.0) + (diff_time.tv_usec/1000.0);
  return time;
}


int appl_timersub(struct timeval* stop_time, struct timeval* start_time, struct timeval* diff_time)
{  
  diff_time->tv_sec  = stop_time->tv_sec  - start_time->tv_sec;
  diff_time->tv_usec = stop_time->tv_usec - start_time->tv_usec;
  return 0;
}









