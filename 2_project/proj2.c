/**
 * @file	proj2.c
 * @brief	2. project IOS  
 * @author	xbarta50
 */

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdbool.h>
#include <pthread.h>
#include <time.h>
#include <errno.h>
#include <semaphore.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/stat.h>
#include <signal.h>
#include <sys/mman.h>

//Global variables
int NO; //input argv[1] for number of oxygen
int NH; //input argv[2] for number of hydrogens
int TI; //input argv[3] for maximal time delay oxygen / hydrogen atom waits before forming molecule
int TB; //input argv[4] for maximal time delay required to create a single molecule
typedef struct {
    int numberofline; //count number of line
    int oxygen; //when forming molecule to check if there is any oxygen to create molecule with
    int hydrogen; //when forming molecule to check if there is 2 hydrogens to create molecule with
    int count; //used at barrier to count if all 3 atoms arrive at barrier
    int molecules; //number of molecule we are creating for print
    int molecules_created; //number of molecule that is created for print
    int create_count; //count if 1 oxygen and 2 hydrogens already printed creating molecule
    int created_count; //count if 1 oxygen and 2 hydrogens already printed molecule created
    int total_Molecules; //number of how many molecules can be created
    int atom_queue; //counter for how many atoms have already printed going to queue
    int count_oxyqueue; //counter for printing O not enough
    int count_hydroqueue; //counter for printing H not enough

    sem_t not_enough; //semphore where molecules are waiting to print not enough after last molecule is creating
    sem_t output; //semaphore for printing out
    sem_t mutex; //semaphore for protecting oxygen and hydrogen counters
    sem_t mutex1; //semaphore for procecting reusable barrier counters
    sem_t hydroQueue; //semaphore where hydrogens are waiting until there is 1 hydrogen and 1 oxygen to create molecule with
    sem_t oxyQueue; //semaphore where oxygens are waiting until there are 2 hydrogens to create molecule with
    sem_t turnstile; //semaphore used in reusable barrier to wait until all 3 atoms arrived
    sem_t turnstile2; //semaphore used in reusable barrier to wait until all 3 atoms arrived
} sharedMem;
sharedMem *shmem; //shared memory segment pointed
FILE *file; //file pointer
// forward declarations
void check_args(int argc, char *argv[]);
void load_args(char* argv[]);
void init();
void Oxygen();
void Hydrogen();
void Barrierwait();
// Clean up functions
void clean();

int main(int argc, char* argv[]) {
    //Functions for the correct number of parameters and if they are only numbers
    check_args(argc, argv);
    //Function that load parameters and checks if they are in a given range 
    load_args(argv);
    //Initialize semaphores, memory and file we will write into
    init();
    //calculates how many molecules will be created
    for (int i = NO, j = NH; i > 0 && j > 1; i--, j-=2){
        shmem->total_Molecules++;
    }
    //Oxygen process
    for(int i=1;i<=NO;++i){
        pid_t pid = fork();
        if(pid < 0){ //forking error occurred
            fprintf(stderr,"Oxygen process creating failed\n");
            clean();
            exit(1);
        }
        if (pid == 0){//Oxygen child process enters its function, at the end of which it kills itself
            Oxygen(i);
        }
    }
    //Hydrogen process
    for(int i=1;i<=NH;++i){
        pid_t pid = fork();
        if(pid < 0){ //forking error occurred
            fprintf(stderr,"Hydrogen process creating failed\n");
            clean();
            exit(1);
        }
        if (pid == 0){//Hydrogen child process enters its function, at the end of which it kills itself
            Hydrogen(i);
        }
    }

    // wait for both generators to exit
    while (wait(NULL) > 0);
    //Release and clear the shared memory
    clean();
    fclose(file);
    return 0;
}
/**
 *  @brief Check if there is 5 arguments and they are only numbers
 */
void check_args(int argc, char* argv[]){
    // check whether 5 arguments have been entered
    if (argc != 5){
        fprintf(stderr, "Wrong number of arguments\n");
        exit(1);
    }
    //check if in num is string
    long ifnumber;
    char *letters;
    for (int i = 1; i < 5; i++) { //check all arguments
        ifnumber = strtol (argv[i], &letters, 10);  //get value of arguments
        if ((*letters != '\0') || (letters==argv[i]) || (ifnumber%1 != 0)){ //check for characters and if its number
			fprintf(stderr, "Arguments are only numbers\n");
			exit(1);
		}
    }   
}
/**
 *  @brief Loads and stores arguments from standard input into pre-prepared global variables
 */
void load_args(char* argv[]){
    //save argumets
    NO = atoi(argv[1]);
    NH = atoi(argv[2]);
    TI = atoi(argv[3]);
    TB = atoi(argv[4]);

    if( NO<=0 || *argv[1] == '\0'){
        fprintf(stderr, "Oxygen must be greater than 0\n");
        exit(1);
    }

    if( NH<=0 || *argv[2] == '\0'){
        fprintf(stderr, "Hydrogen must be greater than 0\n");
        exit(1);
    }

    if( TI<0 || 1000<TI || *argv[3] == '\0'){
        fprintf(stderr, "Third parametr can be 0<=TI<=1000\n");
        exit(1);
    }
    if( TB<0 || 1000<TB || *argv[4] == '\0'){
        fprintf(stderr, "Fourth parametr can be 0<=TB<=1000\n");
        exit(1);
    }
}
/**
 * @brief Function opens file, maps and initialize semaphores
 * @param shmem pointer to the shared memory structure 
 * @param file file where we write
 */
void init(){
    file = fopen("proj2.out", "w");
	if (file == NULL) {
		fprintf(stderr, "Fail to open proj2.out.\n");
		exit(1);
	}

    shmem = mmap(NULL, sizeof(sharedMem), PROT_READ | PROT_WRITE, MAP_SHARED | MAP_ANONYMOUS, -1, 0);
    if (shmem == MAP_FAILED){
        fprintf(stderr, "Fail to initializes shared memory structure\n");
        exit(1);
    }

    shmem->numberofline = 1;
    shmem->molecules = 1;
    shmem->molecules_created = 1;

    if (sem_init(&shmem->output, 1, 1) == -1     || 
        sem_init(&shmem->mutex, 1, 1) == -1      ||
        sem_init(&shmem->mutex1, 1, 1) == -1     ||
        sem_init(&shmem->oxyQueue, 1, 0) == -1   ||
        sem_init(&shmem->hydroQueue, 1, 0) == -1 ||
        sem_init(&shmem->turnstile, 1, 0) == -1  ||
        sem_init(&shmem->turnstile2, 1, 1) == -1 ||
        sem_init(&shmem->not_enough, 1, 0) == -1 ){
        fprintf(stderr, "Fail to initializes the unnamed semaphores\n");
        clean();
        exit(1);
    }
         
}
/**
 * @brief Release memory and destroy semaphores
 */
void clean(){
    if (sem_destroy(&shmem->output) == -1     || 
        sem_destroy(&shmem->mutex) == -1      ||
        sem_destroy(&shmem->mutex1) == -1     ||
        sem_destroy(&shmem->oxyQueue) == -1   ||
        sem_destroy(&shmem->hydroQueue) == -1 ||
        sem_destroy(&shmem->turnstile) == -1  ||
        sem_destroy(&shmem->turnstile2) == -1 ||
        sem_destroy(&shmem->not_enough) == -1 ){
        fprintf(stderr, "Fail to destroy the unnamed semaphores\n");
        exit(1);
    }
    if (munmap(shmem, sizeof(sharedMem))){
        fprintf(stderr, "munmap failed\n");
    }
}
/**
 *  @brief Oxygen function
 *  @param shmem: pointer to the shared memory structure
 */
void Oxygen(int i){
    srand(time(NULL) * getpid()); //randomize the seed

    sem_wait(&shmem->output);
    fprintf(file,"%d: O %d: started\n",shmem->numberofline++,i);
    fflush(file);
    sem_post(&shmem->output);

    usleep(1000 * (rand() % (TI + 1))); //waits for the formation of molecules

    sem_wait(&shmem->output);
    fprintf(file,"%d: O %d: going to queue\n",shmem->numberofline++,i);
    fflush(file);
    shmem->atom_queue++;
    sem_post(&shmem->output);
    if(shmem->atom_queue == NH + NO){ //if you are last atom
        for(int i=0;i<NO - shmem->total_Molecules+NH - shmem->total_Molecules*2;i++){ //open not_enough print as many times as needed
            sem_post(&shmem->not_enough);
        }
        if(shmem->total_Molecules == 0){ //if no molecule will be created
            for(int i = 0;i<NO - shmem->total_Molecules;i++){ //open semaphore for oxyQueue because i open semaphore for oxyQueue after last molecule is creating and there will not be any molecule created
                sem_post(&shmem->oxyQueue); 
            }
        }
        if(shmem->total_Molecules == 0){ //if no molecule will be created
            for(int i = 0;i<NH - shmem->total_Molecules*2;i++){ //open semaphore for hydroQueue because i open semaphore for oxyQueue after last molecule is creating and there will not be any molecule created
                sem_post(&shmem->hydroQueue);
            }
        }
    }

    sem_wait(&shmem->mutex);
    shmem->oxygen += 1;
    if(shmem->hydrogen>=2){ //if there is already at least one oxygen and two hydrogens
        sem_post(&shmem->hydroQueue); //then release them
        sem_post(&shmem->hydroQueue); //then release them
        shmem->hydrogen -= 2;
        sem_post(&shmem->oxyQueue); //then release them
        shmem->oxygen -= 1; 
    } 
    else {
    sem_post(&shmem->mutex);
    }
    sem_wait(&shmem->oxyQueue); //here they wait until 2 hydrogens come
    shmem->count_oxyqueue++;
    //bond  
    if(shmem->count_oxyqueue > shmem->total_Molecules){ //if you are oxygen a you dont have 2 hydrogens to pair with for creating molecule
        sem_wait(&shmem->not_enough); //here wait until last molecule is creating
        sem_wait(&shmem->output);
        fprintf(file,"%d: O %d: not enough H\n",shmem->numberofline++,i);
        fflush(file);
        sem_post(&shmem->output);
        
        exit(0);
    }

    sem_wait(&shmem->output);
    fprintf(file,"%d: O %d: creating molecule %d\n",shmem->numberofline++,i,shmem->molecules);
    fflush(file);
    shmem->create_count++;
    if (shmem->create_count == 3){ //if all 2 hydrogens and 1 oxygen printed creating molecule then increment molecule for next print
        shmem->molecules += 1;
        shmem->create_count = 0;
    }
    if(shmem->molecules > shmem->total_Molecules){ //if you are last molecule that will be created
        if(NO - shmem->total_Molecules){ //if there are any oxygens that dont have enough hydrogens to create molecule with
            for(int i = 0;i<NO - shmem->total_Molecules;i++){ //then release them so they can print not enough
                sem_post(&shmem->oxyQueue);
            }
        }
        if(NH - shmem->total_Molecules*2){ //if there are any hydrogens that dont have oxygens to create molecule with
            for(int i = 0;i<NH - shmem->total_Molecules*2;i++){ //then release them so they can print not enough
                sem_post(&shmem->hydroQueue);
            }
        }
    }
    sem_post(&shmem->output);

    usleep(1000 * (rand() % (TB + 1))); //waits for the formation of molecules
    
    Barrierwait();
    //created


    sem_wait(&shmem->output);
    fprintf(file,"%d: O %d: molecule %d created\n",shmem->numberofline++,i,shmem->molecules_created);
    fflush(file);
    shmem->created_count++;
    if (shmem->created_count == 3){ //if all 2 hydrogens and 1 oxygen printed molecule created then increment molecules_created for next print 
        shmem->molecules_created += 1;
        shmem->created_count = 0;
    }
    sem_post(&shmem->output);

    sem_post(&shmem->mutex);

    exit(0);
}
/**
 *  @brief Hydrogen function
 *  @param shmem: pointer to the shared memory structure
 */
void Hydrogen(int i){
    srand(time(NULL) * getpid()); //randomize the seed

    sem_wait(&shmem->output);
    fprintf(file,"%d: H %d: started\n",shmem->numberofline++,i);
    fflush(file);
    sem_post(&shmem->output);

    usleep(1000 * (rand() % (TI + 1))); //waits for the formation of molecules

    sem_wait(&shmem->output);
    fprintf(file,"%d: H %d: going to queue\n",shmem->numberofline++,i);
    fflush(file);
    shmem->atom_queue++;
    sem_post(&shmem->output);
    if(shmem->atom_queue == NH + NO){ //if you are last atom
        for(int i=0;i<(NO - shmem->total_Molecules)+NH - shmem->total_Molecules*2;i++){ //open not_enough print as many times as needed
            sem_post(&shmem->not_enough);
        }
        if(shmem->total_Molecules == 0){ //if no molecule will be created 
            for(int i = 0;i<NH - shmem->total_Molecules*2;i++){ //open semaphore for hydroQueue because i open semaphore for hydroQueue after molecule is creating and there will not be any molecule created
                sem_post(&shmem->hydroQueue);
            }
        }
        if(shmem->total_Molecules == 0){ //if no molecule will be created
            for(int i = 0;i<NO - shmem->total_Molecules;i++){ //open semaphore for oxyQueue because i open semaphore for oxyQueue after molecule is creating and there will not be any molecule created
                sem_post(&shmem->oxyQueue);
            }
        }
    }

    sem_wait(&shmem->mutex);
    shmem->hydrogen += 1;
    if(shmem->hydrogen>=2 && shmem->oxygen>=1){ //if there is already at least one oxygen and two hydrogens
        sem_post(&shmem->hydroQueue); //then release them
        sem_post(&shmem->hydroQueue); //then release them
        shmem->hydrogen -= 2;
        sem_post(&shmem->oxyQueue); //then release them
        shmem->oxygen -= 1;
    } 
    else {
    sem_post(&shmem->mutex);
    }
    sem_wait(&shmem->hydroQueue); //here they wait until 1 oxygen comes for 2 hydrogens
    shmem->count_hydroqueue++;
    //bond  
    if(shmem->count_hydroqueue > shmem->total_Molecules*2){ //if you are hydrogen a you dont have 1 hydrogen and 1 oxygen to pair with for creating molecule
        sem_wait(&shmem->not_enough); //here wait until last molecule is creating
        sem_wait(&shmem->output);
        fprintf(file,"%d: H %d: not enough O or H\n",shmem->numberofline++,i);
        fflush(file);
        sem_post(&shmem->output);

        
        exit(0);
    }

    //bond()
    sem_wait(&shmem->output);
    fprintf(file,"%d: H %d: creating molecule %d\n",shmem->numberofline++,i,shmem->molecules);
    fflush(file);
    shmem->create_count++;
    if (shmem->create_count == 3){ //if all 2 hydrogens and 1 oxygen printed creating molecule then increment molecule for next print
        shmem->molecules += 1;
        shmem->create_count = 0;
    }
    if(shmem->molecules > shmem->total_Molecules){ //if you are last molecule that will be created
        if(NO - shmem->total_Molecules){ //if there are any oxygens that dont have enough hydrogens to create molecule with
            for(int i = 0;i<NO - shmem->total_Molecules;i++){ //then release them so they can print not enough
                sem_post(&shmem->oxyQueue);
            }
        }
        if(NH - shmem->total_Molecules*2){ //if there are any hydrogens that dont have oxygens to create molecule with
            for(int i = 0;i<NH - shmem->total_Molecules*2;i++){ //then release them so they can print not enough
                sem_post(&shmem->hydroQueue);
            }
        }
    }
    sem_post(&shmem->output);
    
    Barrierwait();

    sem_wait(&shmem->output);
    fprintf(file,"%d: H %d: molecule %d created\n",shmem->numberofline++,i,shmem->molecules_created);
    fflush(file);
    shmem->created_count++;
    if (shmem->created_count == 3){ //if all 2 hydrogens and 1 oxygen printed molecule created then increment molecules_created for next print 
        shmem->molecules_created += 1;
        shmem->created_count = 0;
    }
    sem_post(&shmem->output);

    exit(0);
}

void Barrierwait(){
    sem_wait(&shmem->mutex1);
    shmem->count += 1; 
    if(shmem->count == 3){ //wait for all 3 molecules
        sem_wait(&shmem->turnstile2); 
        sem_post(&shmem->turnstile); 
    }
    sem_post(&shmem->mutex1); 
    sem_wait(&shmem->turnstile); //here molecules wait
    sem_post(&shmem->turnstile); //if you are finally 3 molecules open semaphore for each other
    //critical point

    sem_wait(&shmem->mutex1);
    shmem->count -= 1;
    if(shmem->count == 0){ //wait for all 3 molecules
        sem_wait(&shmem->turnstile);
        sem_post(&shmem->turnstile2);
    }
    sem_post(&shmem->mutex1);
    sem_wait(&shmem->turnstile2); //here molecules wait
    sem_post(&shmem->turnstile2); //if they are finally 3 molecules open semaphore for each other
}