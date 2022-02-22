/*
							FIRST ASSIGNMENT OF REAL TIME OPERATING SYSTEMS

MATTEO CARLONE , MAT. S4652067

First and foremost, I approached the Professor request on paper as in the class during lessons.
We have four tasks J1, J2, J3 and J4,  each of them has different periods respectively:

			J1 80 milliseconds
			J2 100 milliseconds
			J3 160 milliseconds
			J4 200 milliseconds

Already ordered and named according to Rate-Monotonic priorities: 
the lower the period the higher the priority of the Task.
The Tasks can access different resources.

			J1 write on two variables T1T2 and T1T4
			J2 read the variable T1T2 and write on the variable T2T3
			J3 read the variable T2T3 
			J4 read the variable T1T4

Therefore for each Task, it is possible to identify critical sections ( zij ) in which access 
to the memory cell containing T1T2, T1T4 and T2T3 occur.

			J1 = { z11 , z12 }
			J2 = { z21 , z22 }
			J3 = { z31 }
			J4 = { z41 }

I created 3 Semaphores S1, S2 and S3 that will protect those critical zones according to the Priority Ceiling protocol,
that prevents from Deadlocks and Blocking in Cascade by assigning a priority ( c(Si)) corresponding to the highest 
among all tasks that can use the semaphore.

			S1 protects z11 and z21 		| 		c(S1) = P1 , priority of J1        
			S2 protects z22 and z31    		| 		c(S2) = P2 , priority of J2
			S3 protects z21 and z41      	| 		c(S3) = P1 , priority of J1

Before starting scheduling I checked whether the sufficient condition U<Ulub, calculated considering shared resources
with priority ceiling:

Where: 
	- U = Utilization Factor 
	- Ulub = Utilization Factor Lower Upper bound
		
∀i, 1< i < n,   ∑ (from k=1 to i) (WCET[k]/Period[k] + B*i/Period[i]) < i*(2^(1/i) -1)

It was, therefore, necessary to calculate:

	- the WCET (worst-case execution time) by running many times all the tasks in Standalone 
	  Modality and saving the largest run-time, of course for each task, 

	- the B*i, the maximum Blocking time of Task i, The higher block created by a critical section
		in the subset ßi*: 

			- ßij are the critical section of Task j that can block Task i .

	        - ß*i is the Set of Critical Section that blocks Task i .

					ß12 = {z21}, ß13 = {∅}, ß14 = {z41}		|	ß1* = {z21 , z41}
															|
					ß23 = {z31}, ß24 = {z41}				|	ß2* = {z31 , z41}
															|
					ß34 = {z41}								|	ß3* = {z41}
															|	
															|	ß4* = {∅}

I implemented a function called waste_time() to increase the time in which the tasks access 
to the critical sections protected by the semaphores.

*/

// To compile:

// g++ -lpthread <sourcename> -o <executablename>

////Including libraries////

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>
#include <math.h>
#include <sys/types.h>
#include <sys/types.h>

///////////////////////////

//// Colors Definition ////

#define RESET   "\033[0m"
#define BHRED "\e[1;91m"
#define BHGRN "\e[1;92m"
#define BHYEL "\e[1;93m"
#define BHBLU "\e[1;94m"

///////////////////////////

//// Functions Definition ////

void waste_time();

//////////////////////////////

//code of periodic tasks

void task1_code();
void task2_code();
void task3_code();
void task4_code();

//characteristic function of the thread, only for timing and synchronization

//periodic tasks

void *task1( void *);
void *task2( void *);
void *task3( void *);
void *task4( void *);

// initialization of mutexes and conditions

pthread_mutex_t T1T2_Protection; 
pthread_mutex_t T1T4_Protection;
pthread_mutex_t T2T3_Protection;

//// LOOPS Definition ////

#define INNERLOOP 300
#define OUTERLOOP 1000

//////////////////////////

//// Initializing Global variables and arrays to compute the scheduling ////

#define NPERIODICTASKS 4
#define NTASKS NPERIODICTASKS 

long int periods[NTASKS];	// PERIODS ARRAY
struct timespec next_arrival_time[NTASKS];	// NEXT ARRIVAL TIME ARRAY
double WCET[NTASKS];	// WORST CASE EXECUTION TIME ARRAY

pthread_attr_t attributes[NTASKS];
pthread_t thread_id[NTASKS];
struct sched_param parameters[NTASKS];
int missed_deadlines[NTASKS];

int T1T2 , T1T4 , T2T3; 	// TiTj: The variable written by the i' task and read by the j' task 

// Initializing the times of the crtical zones Zij: i is the task of the critical section ,
//                                                  j is the number of critical section

double  time_z11, time_z12, time_z21, time_z22, time_z31, time_z41; 

double Bs1[2],	Bs2[2] , Bs3[1] , Bs4[1] ;

int main()
{
    pthread_mutexattr_t mymutexattr;
	pthread_mutexattr_init(&mymutexattr);
	pthread_mutexattr_setprotocol(&mymutexattr,PTHREAD_PRIO_PROTECT);

  	// set task periods in nanoseconds
	// the first task has period 80 millisecond
	// the second task has period 100 millisecond
	// the third task has period 160 millisecond
	// the fourth task has period 200 millisecond
	
	// Already ordered according to their priority 
	// In RateMonotonic Scheduling the lower the period the higher the priority
	
  	periods[0]= 80000000;  // in nanoseconds
  	periods[1]= 100000000; // in nanoseconds
  	periods[2]= 160000000; // in nanoseconds
  	periods[3]= 200000000; // in nanoseconds

  	struct sched_param priomax; // MINIMUM PRIORITY IN THE SYSTEM

  	priomax.sched_priority=sched_get_priority_max(SCHED_FIFO);

  	struct sched_param priomin; // MAXIMUM PRIORITY IN THE SYSTEM 

  	priomin.sched_priority=sched_get_priority_min(SCHED_FIFO);

	// To set the maximum priority to the current thread you are required to be superuser. 

  	// Using getuid() function I Check that the main thread is executed with superuser privileges
	// before doing anything else.

  	if (getuid() == 0)
    		pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomax);

  	// execute all tasks in standalone modality in order to measure execution times.
  	// Use the computed values to update the worst case execution time of each task.

	for(int i = 0; i< NTASKS; i++){ WCET[i] = 0.0; }
	double WCET_TMP[NTASKS];

	printf(BHBLU "\n    SCHEDULABILITY CHECK \n\n" RESET);

	// executing each task more than once to compute the WCET (Worst Execution Time) in a more accurate way

	for (int n=0; n<100 ; n++){

	  	for(int i = 0 ; i < NTASKS ; i++){

			// initialize time_1 and time_2 required to read the clock
			struct timespec time_1, time_2;
			clock_gettime(CLOCK_REALTIME, &time_1);

	 	     	if (i==0)
				task1_code();
	      		if (i==1)
				task2_code();
	      		if (i==2)
				task3_code();
				if (i==3)
				task4_code();

			clock_gettime(CLOCK_REALTIME, &time_2);

	      	WCET_TMP[i]= 1000000000*(time_2.tv_sec - time_1.tv_sec)+(time_2.tv_nsec-time_1.tv_nsec);

	      	if(WCET_TMP[i]>WCET[i]){

	      		WCET[i] = WCET_TMP[i];

	      	}
	      			
	    }
	  
	}

	//// Initializing U(i) and Ulub(i) ////

    double U1,U2,U3,U4 ;
    double Ulub1,Ulub2,Ulub3,Ulub4;

    ///////////////////////////////////////

    //// Initializing ß*i The Set of Critical Section that blocks Task i ////

    Bs1[0] = time_z21;
    Bs1[1] = time_z41;
    Bs2[0] = time_z31;
   	Bs2[1] = time_z41;
   	Bs3[0] = time_z41;
   	Bs4[0] = 0;

   	//////////////////////////////////////////////////////////////////////////

   	//// Checking the sufficent condition for schedulability ////

    for (int i = 0 ; i < NTASKS ; i++){  // for loop over the tasks

    	printf("\n Worst Case Execution Time %d=%f ", i, WCET[i]);
    	fflush(stdout);

    	switch(i){ 

    		case 0: // Task_1

    			// Compute U1 using the maximum blocking time

    			if(Bs1[0]>Bs1[1]){ 
    				U1 = WCET[0]/periods[0] + Bs1[0]/periods[0];
    			}
    			else{ U1 = WCET[0]/periods[0] + Bs1[1]/periods[0]; }
    			
    			Ulub1 = (i+1)*(pow(2.0,1.0/(i+1))-1);

    			if(U1 > Ulub1){
    				//check the sufficient conditions: if they are not satisfied, exit ;
    				printf("\n U1=%lf     Ulub1=%lf" BHRED " Sufficient Non schedulable Task Set\n" RESET, U1, Ulub1);
    				fflush(stdout);
      				return(-1);
    			}

    		break;

    		case 1: // Task_2

    			// Compute U2 using the maximum blocking time

    			if(Bs2[0]>Bs2[1]){
    				U2 = WCET[0]/periods[0] + WCET[1]/periods[1] + Bs2[0]/periods[1];
    			}
    			else{ U2 = WCET[0]/periods[0] + WCET[1]/periods[1] + Bs2[1]/periods[1]; }
    			
    			Ulub2 = (i+1)*(pow(2.0,1.0/(i+1))-1);
    			if(U2 > Ulub2){
    				//check the sufficient conditions: if they are not satisfied, exit ;
    				printf("\n U2=%lf     Ulub2=%lf" BHRED " Sufficient Condition Not Verified\n" RESET, U2, Ulub2);
    				fflush(stdout);
      				return(-1);
    			}

    		break;

    		case 2: // Task_3

    			U3 = WCET[0]/periods[0] + WCET[1]/periods[1] + WCET[2]/periods[2] + Bs3[0]/periods[2];
    			Ulub3 = (i+1)*(pow(2.0,1.0/(i+1))-1);
    			if(U3 > Ulub3){
    				//check the sufficient conditions: if they are not satisfied, exit ;
    				printf("\n  U3=%lf     Ulub3=%lf" BHRED " Sufficient Condition Not Verified\n" RESET, U3, Ulub3);
    				fflush(stdout);
      				return(-1);
    			}

    		break;

    		case 3: // Task_4

    			U4 = WCET[0]/periods[0] + WCET[1]/periods[1] + WCET[2]/periods[2] + WCET[3]/periods[3] + Bs4[0]/periods[3];
    			Ulub4 = (i+1)*(pow(2.0,1.0/(i+1))-1);
    			if(U4 > Ulub4){
    				//check the sufficient conditions: if they are not satisfied, exit ;
    				printf("\n  U4=%lf     Ulub4=%lf" BHRED " Sufficient Condition Not Verified\n" RESET, U4, Ulub4);
    				fflush(stdout);
      				return(-1);
    			}

    		break;
    	}
    }
  	printf("\n\n" BHGRN "                 Scheduable Task Set\n" RESET);
  	printf("\n   U1=%lf          Ulub1=%lf", U1, Ulub1);
  	printf("\n   U2=%lf          Ulub2=%lf", U2, Ulub2);
  	printf("\n   U3=%lf          Ulub3=%lf", U3, Ulub3);
  	printf("\n   U4=%lf          Ulub4=%lf\n", U4, Ulub4);
  	sleep(2);

  	/////////////////////////////////////////////////////////////

  	// set the minimum priority to the current thread: this is now required because 
	// we will assign higher priorities to periodic threads to be soon created
	// pthread_setschedparam

  	if (getuid() == 0)
    		pthread_setschedparam(pthread_self(),SCHED_FIFO,&priomin);

  
  	// set the attributes of each task, including scheduling policy and priority
  	for (int i =0; i < NPERIODICTASKS; i++)
    	{
		// initialize the attribute structure of task i
      		pthread_attr_init(&(attributes[i]));

		// set the attributes to tell the kernel that the priorities and policies are explicitly chosen,
		// not inherited from the main thread (pthread_attr_setinheritsched) 
      		pthread_attr_setinheritsched(&(attributes[i]), PTHREAD_EXPLICIT_SCHED);
      
		// set the attributes to set the SCHED_FIFO policy (pthread_attr_setschedpolicy)
			pthread_attr_setschedpolicy(&(attributes[i]), SCHED_FIFO);

		// properly set the parameters to assign the priority inversely proportional 
		// to the period
      		parameters[i].sched_priority = priomin.sched_priority+NTASKS - i;

		// set the attributes and the parameters of the current thread (pthread_attr_setschedparam)
      		pthread_attr_setschedparam(&(attributes[i]), &(parameters[i]));
    	}

    //// Setting the Priority Ceiling to the Mutexes ////
    // Where parameters[].sched_priority is ceiling of the semaphore 
    pthread_mutexattr_setprioceiling(&mymutexattr,parameters[0].sched_priority);
    pthread_mutex_init(&T1T2_Protection,&mymutexattr);
    pthread_mutexattr_setprioceiling(&mymutexattr,parameters[0].sched_priority);
    pthread_mutex_init(&T1T4_Protection,&mymutexattr);
    pthread_mutexattr_setprioceiling(&mymutexattr,parameters[1].sched_priority);
    pthread_mutex_init(&T2T3_Protection,&mymutexattr);
    /////////////////////////////////////////////////////

	//delare the variable to contain the return values of pthread_create	
  	int iret[NTASKS];

	//declare variables to read the current time
	struct timespec time_1;
	clock_gettime(CLOCK_REALTIME, &time_1);

  	// set the next arrival time for each task. This is not the beginning of the first
	// period, but the end of the first period and beginning of the next one. 
  	for (int i = 0; i < NPERIODICTASKS; i++)
    	{
		long int next_arrival_nanoseconds = time_1.tv_nsec + periods[i];
		//then we compute the end of the first period and beginning of the next one
		next_arrival_time[i].tv_nsec= next_arrival_nanoseconds%1000000000;
		next_arrival_time[i].tv_sec= time_1.tv_sec + next_arrival_nanoseconds/1000000000;
       		missed_deadlines[i] = 0;
    	}

	// create all threads(pthread_create)
  	iret[0] = pthread_create( &(thread_id[0]), &(attributes[0]), task1, NULL);
  	iret[1] = pthread_create( &(thread_id[1]), &(attributes[1]), task2, NULL);
  	iret[2] = pthread_create( &(thread_id[2]), &(attributes[2]), task3, NULL);
  	iret[4] = pthread_create( &(thread_id[3]), &(attributes[3]), task4, NULL);

  	// join all threads (pthread_join)
  	pthread_join( thread_id[0], NULL);
  	pthread_join( thread_id[1], NULL);
  	pthread_join( thread_id[2], NULL);
  	pthread_join( thread_id[3], NULL);

  	for (int i = 0; i < NTASKS; i++)
    	{
      		printf (" Missed Deadlines Task %d=%d\n", i, missed_deadlines[i]);
		fflush(stdout);
    	}

    pthread_mutexattr_destroy(&mymutexattr);

  	exit(0);
}


// application specific task_1 code
void task1_code()
{
	//print the id of the current task
  	printf(BHRED "1[ " RESET); fflush(stdout);

 	struct timespec time_1, time_2;

 	// Protecting the access to the shared resources using the mutex, lock and unlock
	pthread_mutex_lock(&T1T2_Protection);
	// Getting the duration of of the critical section
	clock_gettime(CLOCK_REALTIME, &time_1);
	// function to waste time to increase the access to the shared resource
	waste_time();
	T1T2 = rand()*rand();
	clock_gettime(CLOCK_REALTIME, &time_2);
	pthread_mutex_unlock(&T1T2_Protection);
	
	// Calculate the effective duration of the critical section
  	time_z11 = 1000000000*(time_1.tv_sec-time_2.tv_sec)+(time_1.tv_nsec-time_2.tv_nsec);

  	struct timespec time_3 , time_4;
		
	pthread_mutex_lock(&T1T4_Protection);
	clock_gettime(CLOCK_REALTIME, &time_3);
	waste_time();
	T1T4 = rand()*rand();
	clock_gettime(CLOCK_REALTIME, &time_4);
	pthread_mutex_unlock(&T1T4_Protection);
		
	time_z12 = 1000000000*(time_3.tv_sec-time_4.tv_sec)+(time_3.tv_nsec-time_4.tv_nsec);
  
  	//print the id of the current task
  	printf(BHRED "]1 " RESET); fflush(stdout);
  	
}

//thread code for task_1 (used only for temporization)
void *task1( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

   	//execute the task one hundred times... it should be an infinite loop (too dangerous)
  	int i=0;
  	for (i=0; i < 100; i++)
    	{
      		// execute application specific code
			task1_code();

			//// checking if we missed deadlines ////
			struct timespec currenttime;
			clock_gettime(CLOCK_REALTIME, &currenttime);
			double currenttimed = currenttime.tv_sec*1000000000 + currenttime.tv_nsec;
			double arrivaltime = next_arrival_time[0].tv_sec*1000000000 + next_arrival_time[0].tv_nsec;
			if(currenttimed>arrivaltime) missed_deadlines[0]++;
			//////////////////////////////////////////

			// sleep until the end of the current period (which is also the start of the
			// new one
			clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[0], NULL);

			// the thread is ready and can compute the end of the current period for
			// the next iteration
	 		
			long int next_arrival_nanoseconds = next_arrival_time[0].tv_nsec + periods[0];
			next_arrival_time[0].tv_nsec= next_arrival_nanoseconds%1000000000;
			next_arrival_time[0].tv_sec= next_arrival_time[0].tv_sec + next_arrival_nanoseconds/1000000000;
    	}
}

void task2_code()
{
	//print the id of the current task
  	printf(BHYEL "2[ " RESET); fflush(stdout);
  	int value_two;
	struct timespec time_1 ,time_2;
	
	pthread_mutex_lock(&T1T2_Protection);
	clock_gettime(CLOCK_REALTIME, &time_1);
	waste_time();
	value_two = T1T2;
	clock_gettime(CLOCK_REALTIME, &time_2);
	pthread_mutex_unlock(&T1T2_Protection);

    time_z21 = 1000000000*(time_1.tv_sec-time_2.tv_sec)+(time_1.tv_nsec-time_2.tv_nsec);

	//print the id of the current task
    
  	struct timespec time_3 , time_4;
  	
	pthread_mutex_lock(&T2T3_Protection);
	clock_gettime(CLOCK_REALTIME, &time_3);
	waste_time();
	T2T3 = rand()*rand();
	clock_gettime(CLOCK_REALTIME, &time_4);
	pthread_mutex_unlock(&T2T3_Protection);

    time_z22 = 1000000000*(time_3.tv_sec-time_4.tv_sec)+(time_3.tv_nsec-time_4.tv_nsec);

  	printf(BHYEL "]2 " RESET); fflush(stdout);
    
}


void *task2( void *ptr )
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	int i=0;
  	for (i=0; i < 100; i++)
    	{
      		task2_code();

      		struct timespec currenttime;
			clock_gettime(CLOCK_REALTIME, &currenttime);
			double currenttimed = currenttime.tv_sec*1000000000 + currenttime.tv_nsec;
			double arrivaltime = next_arrival_time[1].tv_sec*1000000000 + next_arrival_time[1].tv_nsec;
			if(currenttimed>arrivaltime) missed_deadlines[1]++;

			clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[1], NULL);
			long int next_arrival_nanoseconds = next_arrival_time[1].tv_nsec + periods[1];
			next_arrival_time[1].tv_nsec= next_arrival_nanoseconds%1000000000;
			next_arrival_time[1].tv_sec= next_arrival_time[1].tv_sec + next_arrival_nanoseconds/1000000000;
    	}
}

void task3_code()
{
	//print the id of the current task

  	printf(BHGRN "3[ " RESET); fflush(stdout);
  	int value_three;
	struct timespec time_1, time_2;
	
    pthread_mutex_lock(&T2T3_Protection);
    clock_gettime(CLOCK_REALTIME, &time_1);	
    waste_time();
	value_three = T2T3;  
	clock_gettime(CLOCK_REALTIME, &time_2);
	pthread_mutex_unlock(&T2T3_Protection); 
				
	//print the id of the current task
	time_z31 = 1000000000*(time_1.tv_sec-time_2.tv_sec)+(time_1.tv_nsec-time_2.tv_nsec);
  	printf(BHGRN "]3 " RESET); fflush(stdout);
}

void *task3( void *ptr)
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	int i=0;
  	for (i=0; i < 100; i++)
    	{
      		task3_code();

      		struct timespec currenttime;
			clock_gettime(CLOCK_REALTIME, &currenttime);
			double currenttimed = currenttime.tv_sec*1000000000 + currenttime.tv_nsec;
			double arrivaltime = next_arrival_time[2].tv_sec*1000000000 + next_arrival_time[2].tv_nsec;
			if(currenttimed>arrivaltime) missed_deadlines[2]++;

			clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[2], NULL);
			long int next_arrival_nanoseconds = next_arrival_time[2].tv_nsec + periods[2];
			next_arrival_time[2].tv_nsec= next_arrival_nanoseconds%1000000000;
			next_arrival_time[2].tv_sec= next_arrival_time[2].tv_sec + next_arrival_nanoseconds/1000000000;
    }
}

void task4_code()
{
	
	//print the id of the current task
  	printf(BHBLU "4[ " RESET); fflush(stdout);
  	int value_four;
	struct timespec time_1, time_2;
	pthread_mutex_lock(&T1T4_Protection);
	clock_gettime(CLOCK_REALTIME, &time_1);	
	waste_time();
	value_four = T1T4;
	clock_gettime(CLOCK_REALTIME, &time_2);	
	pthread_mutex_unlock(&T1T4_Protection);
			
	//print the id of the current task
	time_z41 = 1000000000*(time_1.tv_sec-time_2.tv_sec)+(time_1.tv_nsec-time_2.tv_nsec);
  	printf(BHBLU "]4 " RESET); fflush(stdout);

}


void *task4( void *ptr )
{
	// set thread affinity, that is the processor on which threads shall run
	cpu_set_t cset;
	CPU_ZERO (&cset);
	CPU_SET(0, &cset);
	pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cset);

	int i=0;
  	for (i=0; i < 100; i++)
    	{
      		task4_code();

      		struct timespec currenttime;
			clock_gettime(CLOCK_REALTIME, &currenttime);
			double currenttimed = currenttime.tv_sec*1000000000 + currenttime.tv_nsec;
			double arrivaltime = next_arrival_time[3].tv_sec*1000000000 + next_arrival_time[3].tv_nsec;
			if(currenttimed>arrivaltime) missed_deadlines[3]++;
			clock_nanosleep(CLOCK_REALTIME, TIMER_ABSTIME, &next_arrival_time[3], NULL);
			long int next_arrival_nanoseconds = next_arrival_time[3].tv_nsec + periods[3];
			next_arrival_time[3].tv_nsec= next_arrival_nanoseconds%1000000000;
			next_arrival_time[3].tv_sec= next_arrival_time[3].tv_sec + next_arrival_nanoseconds/1000000000;
    	}
}

void waste_time(){

	//This function compute a double loop with random computation only required to waste time
	int i,j;
	double number = 0;;
  	for (i = 0; i < OUTERLOOP; i++)
    	{
      		for (j = 0; j < INNERLOOP; j++)
      		{
			number += 1;
    		}
  		}

}