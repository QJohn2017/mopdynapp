/*

Important:

1. When testing OLD Mopac agianst NEW It is necessary to keep correspondence between *.dynput and *.mopac files 
2. It seems that all quantum atoms must go first in the input - should be ckecked further
modified by H.B. Zhang 2013.11.12 subroutine callmopac_SCF2 and WriteMopacInputM_SCF2 are added

*/

#include "c_unix.h" 
#include "varglobal.h"
int   NAM_LENGTH=80;
char  MopacType[80];
char  MopdynType[80];
char  MopacSolver[80];
char namMopacInput[80];
char namMopacOutput[80];
char namMopacBase[80];
char namMopacAux[80];
char namMinDist[80];
static char namDynMolBase[80];
static char namDynMolDynput[80];
static char namDynMolBonds[80];
static char namDynMolCharges[80];
static char namDynMolHomo[80];
static int   MopacVersion=-1;
static int   MopIterations=0;
static int   TypeofCalculations=-1;
static char  MopCommand[80];
void MinDist(double X[], double Y[], double Z[], int AtomType[],int* Number);
void WriteMopacInput(char* MopacType,double X[], double Y[], double Z[], int AtomType[],int* Number);
void WriteMopacInputM(char* MopacType,double X[], double Y[], double Z[], int AtomType[],int* Number);
void ExecuteMopac(char* namMopacInput);

void ExtractData(char* namMopacOutput, int* Number, double FMOP[], double* ZPOT); //xiaoqing
void ExtractDataAux(char* namMopacAux, int* Number, double FMOP[], double* ZPOT); //xiaoqing
void strrep(char *str, char oldchar, char newchar);

void callmopac(int* nk, double X[], double Y[], double Z[], int AtomType[], double FMOP[], double* ZPOT)
 {  
//  MinDist(X,Y,Z,AtomType,nk);
//  WriteMopacInput(MopacType,X,Y,Z,AtomType,nk);
/* Enable to account for situations with overlaping H-H atoms */
  WriteMopacInputM(MopacType,X,Y,Z,AtomType,nk);
  ExecuteMopac(namMopacInput); 
  //ExtractData(namMopacOutput,nk,FMOP,ZPOT);
  ExtractDataAux(namMopacAux,nk,FMOP,ZPOT);
  
 }
 
 
void callmopac_SCF2(int* nk, double X[], double Y[], double Z[], int AtomType[], double FMOP[], double* ZPOT)
 {  
//  MinDist(X,Y,Z,AtomType,nk);
//  WriteMopacInput(MopacType,X,Y,Z,AtomType,nk);
/* Enable to account for situations with overlaping H-H atoms */
  WriteMopacInputM_SCF2(MopacType,X,Y,Z,AtomType,nk);
  ExecuteMopac(namMopacInput); 
  //ExtractData(namMopacOutput,nk,FMOP,ZPOT);
  ExtractDataAux(namMopacAux,nk,FMOP,ZPOT);
  
 } 
 
void ExtractDataAux(char* namMopacAux, int* Number, double FMOP[], double* ZPOT) //xiaoqing
{
FILE *fp;
char junk[NAM_LENGTH];
double lfjunk, enelec, ennucl;
int i=0,i1=0; 
char delims[] = "=";
char* result = NULL;
i=0;
i1=0;

 sprintf( namMopacAux,"mopac.aux"); 

 if ( ( fp = fopen(namMopacAux, "r") ) == NULL)
  {
    printf("%s Mopac output .aux file could not be opened\n", namMopacOutput );
    exit( 0 );
  }

/* Extract the potential energy - Total energy as defined in MOPAC */
     fscanf( fp, "%s", junk );  
     while ( strncmp( junk, "ENERGY_ELECTRONIC", 16 ) != 0 )
     {  
      fscanf( fp, "%s", junk );      
	 }
	 result = strtok( junk, delims );
	 result = strtok( NULL, delims );
	 strrep(result,'D','E');
	 enelec = atof(result);

     fscanf( fp, "%s", junk );         
     if ( strncmp( junk, "ENERGY_NUCLEAR", 14 ) != 0 )
	 {
       printf("Nuclear energy is not found in MOPAC output .aux file. Exiting ...\n");
       exit(0);
     }
	 result = strtok( junk, delims );
	 result = strtok( NULL, delims );
	 strrep(result,'D','E');
	 ennucl = atof(result);

     *ZPOT = enelec + ennucl;
     //printf("ENERGY_ELECTRONIC IS %e\n\n",enelec);
     //printf("ENERGY_NUCLEAR IS %e\n\n",ennucl);

Begin: 
 
    while ( strncmp( junk, "GRADIENTS:KCAL/MOL/ANGSTROM", 12 ) != 0 )
	{  
      fscanf( fp, "%s", junk );
    }

Middle:

	for (i = 0; i < *Number; i++)
	{
		fscanf( fp, "%lf", &lfjunk); 
		FMOP[i]=lfjunk;

		fscanf( fp, "%lf", &lfjunk); 
		FMOP[i+*Number]=lfjunk;

		fscanf( fp, "%lf", &lfjunk); 
		FMOP[i+(*Number)*2]=lfjunk;

	}
			
Ending:	       
    
    fclose( fp );
}

void strrep(char *str, char oldchar, char newchar) 
{
	char *pos = NULL;
    if (newchar == oldchar) {
        return;
    }
    pos = strchr(str, oldchar);
    while (pos != NULL)  {
        *pos = newchar;
        pos = strchr(pos + 1, oldchar);
    }
}

void ExtractData(char* namMopacOutput, int* Number, double FMOP[], double* ZPOT)
{  
FILE *fp;
char junk[NAM_LENGTH];
float fjunk;
int i=0,i1=0; 
i=0;
i1=0;
 if ( ( fp = fopen(namMopacOutput, "r") ) == NULL)
  {
    printf("%s Mopac output file could not be opened\n", namMopacOutput );
    exit( 0 );
  }

/* Extract the potential energy - Total energy as defined in MOPAC */

  
     fscanf( fp, "%s", junk );  

while ( strcmp( junk, "TOTAL" ) != 0 ) {  
      fscanf( fp, "%s", junk );      
    }
    
     fscanf( fp, "%s", junk );         
     if ( strcmp( junk, "ENERGY" ) != 0 ){
     printf("Total energy is not found in MOPAC output file. Exiting ...\n");
     exit(0);
     }
     fscanf( fp, "%s", junk ); 
     fscanf( fp, "%e", &fjunk ); 
     *ZPOT=fjunk;
//     printf("ENERGY IS %e\n\n",fjunk);

Begin: 
 
while ( strcmp( junk, "FINAL" ) != 0 ) {  
      fscanf( fp, "%s", junk );      
    }
      fscanf( fp, "%s", junk );   
     if ( strcmp( junk, "POINT" ) != 0 )    goto Begin;
       fscanf( fp, "%s", junk );  
       if ( strcmp( junk, "AND" ) != 0 )         goto Begin;
         fscanf( fp, "%s", junk );  
          if ( strcmp( junk, "DERIVATIVES" ) != 0 )    goto Begin;
Middle:	  
              fscanf( fp, "%s", junk );  	     
     		while ( strcmp( junk, "CARTESIAN" ) != 0 ) {  
		       fscanf( fp, "%s", junk );      
		        if( strcmp( junk, "CHEMICAL" ) == 0 ) goto Ending;
	       }	        
	          fscanf( fp, "%s", junk ); 
		     fscanf( fp, "%e", &fjunk ); 
		        fscanf( fp, "%e", &fjunk); 
			FMOP[i]=fjunk;		
			i1++;
			if(i1%3 != 0)
			{
			i=i+*Number;				
			}else{
			i=i-(*Number)*2;
                        i++;
			}			
				
			goto Middle;
			
Ending:	       
    
fclose( fp );
}

void WriteMopacInput(char* MopacType,double X[], double Y[], double Z[], int AtomType[],int* Number)
{
FILE * pFile;
int dummy=1, i=0;
 if ( ( pFile = fopen(namMopacInput, "w") ) == NULL)
  {
    printf("%s Mopac input file could not be opened\n",namMopacInput);
    exit( 0 );
  }
 if(MopIterations == 0){
   fprintf( pFile, "%s %s denout xyz gradients dcart geo-ok\n",MopacType,MopCommand);
 }else{
   fprintf( pFile, "%s %s denout oldens xyz gradients dcart geo-ok\n",MopacType,MopCommand);

 }
 fprintf( pFile, "\n");	
 fprintf( pFile, "\n");	 
  for(i=0;i<*Number;i++){
        if(AtomType[i]==1){
 	 dummy=1;
	}else{
	 dummy=1;
        } 
 fprintf( pFile, "%2d %16.12f %2d %16.12f %2d %16.12f %2d \n", AtomType[i], X[i], dummy, Y[i], dummy, Z[i], dummy); 
 }
 fclose(pFile);
}


void WriteMopacInputM(char* MopacType,double X[], double Y[], double Z[], int AtomType[],int* Number)
{
FILE * pFile;
int dummy=1, i=0, j;
int Match[500];

 for(i=0;i<*Number;i++) Match[i]=0; 
 
 for(i=0;i<*Number;i++)
 if( AtomType[i] == 1){
     
     for(j=i+1;j<*Number;j++) 
      if( (X[i] == X[j]) && (Y[i] == Y[j]) && (Z[i] == Z[j])){
      Match[i]++;
      AtomType[j]=0;      
     }
 }
 
  for(i=0;i<*Number;i++){    
//  printf("1 i=%d,AtomType[i]=%d,Match[i]=%d\n",i,AtomType[i],Match[i]);
  if( AtomType[i] == 1){
    
	  if (Match[i] == 0){
//      Make dummy atom to be oxygen
//	    AtomType[i] = 8;	  
	    AtomType[i] = 1; 
	  }else if (Match[i] == 1){
/* Dummy atom is oxygen */
	   AtomType[i] = 8;     
	  }else if (Match[i] == 2){
/* Dummy atom is nitrogen */
	   AtomType[i] = 7;      
	  }else if (Match[i] == 3){
/* Dummy atom is carbon */
	   AtomType[i] = 6;     
	 }else{
	  AtomType[i]=AtomType[0];          
//	  printf("Valency of dummy atom is %d. This valency is not supported. Exiting \n",(Match[i]+1));
//	  exit(0);
 	}
 } 
//   printf("2 i=%d,AtomType[i]=%d,Match[i]=%d\n",i,AtomType[i],Match[i]);

 }

 if ( ( pFile = fopen(namMopacInput, "w") ) == NULL)
  {
    printf("%s Mopac input file could not be opened\n",namMopacInput);
    exit( 0 );
  }
 if(MopIterations == 0){
   fprintf( pFile, "%s %s denout xyz gradients dcart geo-ok\n",MopacType,MopCommand);
 }else{
   fprintf( pFile, "%s %s denout oldens xyz gradients dcart geo-ok\n",MopacType,MopCommand);

 }
 
 fprintf( pFile, "\n");	
 fprintf( pFile, "\n");	 
 
 
  for(i=0;i<*Number;i++)
     if( AtomType[i] != 0){   
        
  fprintf( pFile, "%2d %16.12f %2d %16.12f %2d %16.12f %2d \n", AtomType[i], X[i], dummy, Y[i], dummy, Z[i], dummy); 
  }
 fclose(pFile);
}

void WriteMopacInputM_SCF2(char* MopacType,double X[], double Y[], double Z[], int AtomType[],int* Number)
{
FILE * pFile;
int dummy=1, i=0, j;
int Match[500];

 for(i=0;i<*Number;i++) Match[i]=0; 
 
 for(i=0;i<*Number;i++)
 if( AtomType[i] == 1){
     
     for(j=i+1;j<*Number;j++) 
      if( (X[i] == X[j]) && (Y[i] == Y[j]) && (Z[i] == Z[j])){
      Match[i]++;
      AtomType[j]=0;      
     }
 }
 
  for(i=0;i<*Number;i++){    
//  printf("1 i=%d,AtomType[i]=%d,Match[i]=%d\n",i,AtomType[i],Match[i]);
  if( AtomType[i] == 1){
    
	  if (Match[i] == 0){
//      Make dummy atom to be oxygen
//	    AtomType[i] = 8;	  
	    AtomType[i] = 1; 
	  }else if (Match[i] == 1){
/* Dummy atom is oxygen */
	   AtomType[i] = 8;     
	  }else if (Match[i] == 2){
/* Dummy atom is nitrogen */
	   AtomType[i] = 7;      
	  }else if (Match[i] == 3){
/* Dummy atom is carbon */
	   AtomType[i] = 6;     
	 }else{
	  AtomType[i]=AtomType[0];          
//	  printf("Valency of dummy atom is %d. This valency is not supported. Exiting \n",(Match[i]+1));
//	  exit(0);
 	}
 } 
//   printf("2 i=%d,AtomType[i]=%d,Match[i]=%d\n",i,AtomType[i],Match[i]);

 }

 if ( ( pFile = fopen(namMopacInput, "w") ) == NULL)
  {
    printf("%s Mopac input file could not be opened\n",namMopacInput);
    exit( 0 );
  }// NOT FINISHED
 if(MopIterations == 0){
   fprintf( pFile, "PM6 uhf bonds 1SCF SCFCRT=1.D-2 AUX(PRECISION=9) denout xyz gradients dcart geo-ok\n");
 }else{
   fprintf( pFile, "PM6 uhf bonds 1SCF SCFCRT=1.D-2 AUX(PRECISION=9) denout oldens xyz gradients dcart geo-ok\n");

 }
 
 fprintf( pFile, "\n");	
 fprintf( pFile, "\n");	 
 
 
  for(i=0;i<*Number;i++)
     if( AtomType[i] != 0){   
        
  fprintf( pFile, "%2d %16.12f %2d %16.12f %2d %16.12f %2d \n", AtomType[i], X[i], dummy, Y[i], dummy, Z[i], dummy); 
  }
 fclose(pFile);
}

void MinDist(double X[], double Y[], double Z[], int AtomType[],int* Number)
{
FILE * pFile;
int i=0, j=0;
double R, Dist[*Number];
if ( ( pFile = fopen(namMinDist, "a") ) == NULL)
  {
    printf("%s MinDist file could not be opened\n",namMinDist);
    exit( 0 );
}

for(i=0;i<*Number;i++){
 Dist[i]=100;
 R=1000;
  for(j=0;j<*Number;j++){
    if (i != j)
     R=sqrt( (X[i]-X[j])*(X[i]-X[j]) + (Y[i]-Y[j])*(Y[i]-Y[j]) + (Z[i]-Z[j])*(Z[i]-Z[j]) );  
     if (Dist[i] > R ) Dist[i]=R;
  }
   
}
  for(i=0;i<*Number;i++)
   fprintf( pFile, "%16.12f \n", Dist[i] );   
     fprintf( pFile, "\n\n", Dist[i] );  
    fclose(pFile);
}


void ExecuteMopac(char* namMopacInput)
{
  char cname[100];
  int istat;

  if( TypeofCalculations >= 50){
  sprintf( cname,"%s %s \n",MopacSolver,namMopacBase);   
  }else{
  sprintf( cname,"wine ./MOPAC_7.2.exe %s \n",namMopacInput); 
  }
  istat =  system(cname);  
  if(istat != 0){
  printf("Mopac run is failed. istat=%d. Exiting... \n",istat);
//  exit(0);
  }
  sprintf( cname,"rm -f mopac.ar* \n"); 
  istat =  system(cname);  
  MopIterations++;
}


void readmopactype(int* imopcall, char Molnam[])
{
FILE * fp;
 char aline[80];
 char junk[NAM_LENGTH],junk1[NAM_LENGTH],junk2[NAM_LENGTH],junk3[NAM_LENGTH],junk4[NAM_LENGTH],junk5[NAM_LENGTH],junk6[NAM_LENGTH]; 
 int i;
   for(i=0;i<NAM_LENGTH;i++){   
     if (Molnam[i] == ' ') break;
       namDynMolBase[i]=Molnam[i];
   }
 
 strcpy(MopCommand,"uhf");
   sprintf( namDynMolDynput,"%s.dynput",namDynMolBase);   
 
 if ( ( fp = fopen(namDynMolDynput, "r") ) == NULL)
  {
    printf("Input file %s could not be opened\n",namDynMolDynput);
    exit( 0 );
  }
    fgets( aline, 80, fp );
    fscanf( fp, "%s %s", junk,junk1);    
    
   if ( strcmp( junk, "MOPAC" ) == 0 ){ 
    
    if ( strcmp( junk1, "VERSION" ) == 0 ){
    *imopcall = 1;
    }else if ( strcmp( junk1, "6" ) == 0 ){ 
    *imopcall = 0;
     MopacVersion=6;
     strcpy(MopCommand,"uhf");
    }else if ( strcmp( junk1, "7" ) == 0 ){ 
    *imopcall = 0;
     MopacVersion=7; 
     strcpy(MopCommand,"C.I.=(5,5) pulay 1SCF");    
    }else{    
     printf("Foramt of Input file %s is not correct - Line 2",namDynMolDynput);
     exit(0);
    }  
   }else  if ( strcmp( junk, "MOPDYN" ) == 0 ){    
    MopacVersion=7; 
      fscanf( fp, "%s", MopdynType); 
       if ( strcmp( MopdynType, "RATES" ) == 0 ){  
       *imopcall = 50;
       }else if  ( strcmp( MopdynType, "MORPHOLOGY" ) == 0 ){  
       *imopcall = 100;
       }else if  ( strcmp( MopdynType, "HEATTRANSFER" ) == 0 ){  
       *imopcall = 500;       
        fscanf( fp, "%s = %s %s = %s %s = %s %s = %s %s = %s %s = %s", junk,junk1,junk,junk2,junk,junk3,junk,junk4,junk,junk5,junk,junk6);
//        printf("junk1=%s, junk2=%s junk3=%s junk4=%s junk5=%s junk6=%s\n",junk1,junk2,junk3,junk4,junk5,junk6);
	if ( strcmp( junk1, "T" ) == 0 ){
             ActiveBndr[0]=1;
	 }else{
             ActiveBndr[0]=0;
	 }
	if ( strcmp( junk2, "T" ) == 0 ){
             ActiveBndr[1]=1;
	 }else{
             ActiveBndr[1]=0;
	 }
	if ( strcmp( junk3, "T" ) == 0 ){
             ActiveBndr[2]=1;
	 }else{
             ActiveBndr[2]=0;
	 }	 
	if ( strcmp( junk4, "T" ) == 0 ){
             ActiveBndr[3]=1;
	 }else{
             ActiveBndr[3]=0;
	 }	
	 if ( strcmp( junk5, "T" ) == 0 ){
             ActiveBndr[4]=1;
	 }else{
             ActiveBndr[4]=0;
	 }
	 if ( strcmp( junk6, "T" ) == 0 ){
             ActiveBndr[5]=1;
	 }else{
             ActiveBndr[5]=0;
	 }
	
printf("ActiveBndr[0]=%d,ActiveBndr[1]=%d,ActiveBndr[2]=%d,ActiveBndr[3]=%d,ActiveBndr[4]=%d,ActiveBndr[5]=%d",ActiveBndr[0],ActiveBndr[1],ActiveBndr[2],ActiveBndr[3],ActiveBndr[4],ActiveBndr[5]);

       }else{
       printf("Foramt of Input file %s is not correct - Line 2 / 3d Word %s",namDynMolDynput,MopdynType);
       exit(0);
       }
         while ( strcmp( junk, "RUN" ) != 0 ) {
          fscanf( fp, "%s", junk );
         }
	  fgets( MopacSolver, 80, fp );
	  for(i=0;i<80;i++)
          if(MopacSolver[i] == '\n'){
          MopacSolver[i] = ' ';
          break;     
     }    
   }else{
     printf("Foramt of Input file %s is not correct - Line 2 / 1st Word",namDynMolDynput);
     exit(0);
   }

/***/

  
  
/***1st Version of DYNPUT FILE***********************/
    if(*imopcall == 0){    
    fscanf( fp, "%s", MopacType);
    printf ("Mopac version is = %d\n",MopacVersion);
    printf ("MopacType = %s\n",MopacType);
    
    if( (MopacVersion != 6 ) && (MopacVersion != 7 ) ){
     printf("Only MOPAC 6 & MOPAC 7 are supported. Exiting ...\n\n");  
     exit(0);  
    }
    
    if( (MopacVersion == 6 ) && (strcmp(MopacType, "PM6" ) == 0 ) ){
    printf("MOPAC 6 does not support PM6 parameterization. Exiting ...\n\n");
    exit(0);
    }   
    
    }else if(*imopcall == 1){
    
/***2nd Version of DYNPUT FILE***********************/      
   
    fscanf( fp, "%d", &MopacVersion);    
    fscanf( fp, "%s", MopacType);
    printf ("Mopac version is = %d\n",MopacVersion);
    printf ("MopacType = %s\n",MopacType);
    
    if( (MopacVersion != 6 ) && (MopacVersion != 7 ) ){
     printf("Only MOPAC 6 & MOPAC 7 are supported. Exiting ...\n\n");  
     exit(0);  
    }
    
    if( (MopacVersion == 6 ) && (strcmp(MopacType, "PM6" ) == 0 ) ){
    printf("MOPAC 6 does not support PM6 parameterization. Exiting ...\n\n");
    exit(0);
    }    
    fgets( aline, 80, fp );  
    fgets( MopCommand, 80, fp );    
    for(i=0;i<80;i++)
     if(MopCommand[i] == '\n'){
     MopCommand[i]=' ';
     break;     
     }        
    }else if(*imopcall == 50){
/***1st Version of DYNPUT FILE FOR COMPUTING RATES ***/   
    while ( strcmp( junk, "PARAMETERIZATION" ) != 0 ) {
          fscanf( fp, "%s", junk );
         } 
     fscanf( fp, "%s", MopacType);
      printf ("Mopac version is = %d\n",MopacVersion);
       printf ("MopacType = %s\n",MopacType);
    
    fgets( aline, 80, fp );  
    fgets( MopCommand, 80, fp );    
    for(i=0;i<80;i++)
     if(MopCommand[i] == '\n'){
     MopCommand[i]=' ';
     break;     
     }    
    }else if(*imopcall == 100){
/***1st Version of DYNPUT FILE FOR INTEGRATING WITH KMC***/    
   while ( strcmp( junk, "PARAMETERIZATION" ) != 0 ) {
          fscanf( fp, "%s", junk );
         } 
     fscanf( fp, "%s", MopacType);
      printf ("Mopac version is = %d\n",MopacVersion);
       printf ("MopacType = %s\n",MopacType);
    
    fgets( aline, 80, fp );  
    fgets( MopCommand, 80, fp );    
    for(i=0;i<80;i++)
     if(MopCommand[i] == '\n'){
     MopCommand[i]=' ';
     break;     
     }    
    }else if(*imopcall == 500){
/***1st Version of DYNPUT FILE FOR MODELING THERMAL TRANSFER***/    
   while ( strcmp( junk, "PARAMETERIZATION" ) != 0 ) {
          fscanf( fp, "%s", junk );
         } 
     fscanf( fp, "%s", MopacType);
      printf ("Mopac version is = %d\n",MopacVersion);
       printf ("MopacType = %s\n",MopacType);
    
    fgets( aline, 80, fp );  
    fgets( MopCommand, 80, fp );    
    for(i=0;i<80;i++)
     if(MopCommand[i] == '\n'){
     MopCommand[i]=' ';
     break;     
     }    
    }
/********************************************************/    
    
    fclose(fp);   
/* Define Mopav input and output files */     
 

   sprintf( namMopacBase,"mopac");  
   sprintf( namMopacInput,"mopac.dat"); 
   sprintf( namMopacOutput,"mopac.out"); 
   sprintf( namMopacAux,"mopac.aux"); 
   sprintf( namDynMolBonds,"%s.bond",namDynMolBase);   
   sprintf( namDynMolCharges,"%s.charge",namDynMolBase);   
   sprintf( namDynMolHomo,"%s.homo",namDynMolBase);
/**************************************/   
//   sprintf( namMinDist,"mindist.dat");  
/**************************************/
TypeofCalculations=*imopcall;

}


void getmopacversion(int* version)
{
*version=MopacVersion;
}


void getbonds(int* BondAtomNum, double* Otime, int* ImageAtomNum)
{
	FILE *fp, *fp1;
	char junk[81],SymbRow[81],aline[81],aline1[81];
	int i=0, NumbRow,NNN;
	
	sprintf( namMopacOutput,"mopac.out");
	NNN = *BondAtomNum;
	
	if ( ( fp = fopen(namMopacOutput, "r") ) == NULL)
	{
		printf("%s Mopac output file could not be opened\n", namMopacOutput );
		exit( 0 );
	}
	
	if ( ( fp1 = fopen(namDynMolBonds, "a") ) == NULL)
	{
		printf("%s Mopdyn output file for bond orders could not be opened\n", namDynMolBonds );
		exit( 0 );
	}
	
	fprintf( fp1, "Simulation time = %e\n",*Otime);     
/******************************************************************************/  
/******** Extract&Print bond orders and valencies as defined in MOPAC *********/
/*****************************************************************************/
/************************ Find "BOND" in mopac.out file ***********************/
	fscanf( fp, "%s", junk );
	while ( strcmp( junk, "BOND" ) != 0 )
	{
		fscanf( fp, "%s", junk );      
    }    
/*********************** Read the rest of the line ***************************/
	fgets( aline, 81, fp ); 
	fprintf( fp1, "%s%s\n",junk,aline); 
/*********************** Print title of the table *****************************/     
	
    while ( !feof(fp) )
	{
		fgets( aline, 81, fp ); 
		i++;
		
		NumbRow =0;
		sscanf( aline, "%s %d", SymbRow, &NumbRow);
		
		if ( strcmp( SymbRow, "TOTAL" ) == 0 )
		{
			goto Ending;
		}
		
		if ( NumbRow > NNN )
		{
			fgets( aline1, 81, fp ); 
			i++;
		}
		else
		{
			if ( aline[0] == '\n' )
			{
//				fprintf( fp1, "%s",aline);
			}
			else
			{
				fprintf( fp1, "%s",aline);
			}
		}
//		printf( "%s  %d  %d\n", SymbRow, NumbRow, i); 
	}

Ending:
	fprintf( fp1,"\n\n\n");
	fclose( fp );
	fclose( fp1 );
} 


void getcharges(int* BondAtomNum, double* Otime, int* ImageAtomNum)
{
	FILE *fp, *fp1;
	char junk[81],SymbRow[81],aline[81];
	int i=0;
	
	sprintf( namMopacOutput,"mopac.out");
	
	if ( ( fp = fopen(namMopacOutput, "r") ) == NULL)
	{
		printf("%s Mopac output file could not be opened\n", namMopacOutput );
		exit( 0 );
	}
	
	if ( ( fp1 = fopen(namDynMolCharges, "a") ) == NULL)
	{
		printf("%s Mopdyn output file for charges could not be opened\n", namDynMolBonds );
		exit( 0 );
	}
	
	fprintf( fp1, "Simulation time = %e\n",*Otime);     
/******************************************************************************/  
/******** Extract&Print charges as defined in MOPAC *********/
/*****************************************************************************/
/************************ Find "BOND" in mopac.out file ***********************/
		fscanf( fp, "%s", junk );
		while ( strcmp( junk, "CHARGES" ) != 0 )
		{
			fscanf( fp, "%s", junk );
		}  
/*********************** Read the rest of the line ***************************/
	fgets( aline, 81, fp ); 
	fprintf( fp1, "NET ATOMIC %s%s",junk,aline); 
/*********************** Print title of the table *****************************/     
	fgets( aline, 81, fp ); 
	//fprintf( fp1, "%s",aline);
	fgets( aline, 81, fp ); 
	fprintf( fp1, "%s",aline);

    while ( !feof(fp) )
	{
		fgets( aline, 81, fp ); 
		i++;
		
		sscanf( aline, "%s", SymbRow);
		
		if ( strcmp( SymbRow, "CARTESIAN" ) == 0 )
		{
			goto Ending;
		}
		
		if ( aline[0] == '\n' )
		{
			goto Ending;
		}
		else
		{
			fprintf( fp1, "%s",aline);
		}
	}

Ending:
	fprintf( fp1,"\n\n\n");
	fclose( fp );
	fclose( fp1 );
} 

void gethomo(double* Otime)
{
	FILE *fp, *fp1;
	char junk[81];
    float fjunk;
	
	sprintf( namMopacOutput,"mopac.out");
	
	if ( ( fp = fopen(namMopacOutput, "r") ) == NULL)
	{
		printf("%s Mopac output file could not be opened\n", namMopacOutput );
		exit( 0 );
	}
	
	if ( ( fp1 = fopen(namDynMolHomo, "a") ) == NULL)
	{
		printf("%s Mopdyn output file for homo could not be opened\n", namDynMolBonds );
		exit( 0 );
	}
	
    fscanf( fp, "%s", junk );  
    while ( strcmp( junk, "LUMO" ) != 0 )
	{  
      fscanf( fp, "%s", junk );      
    }

	fprintf( fp1, "%e",*Otime);

    fscanf( fp, "%s", junk );      
    fscanf( fp, "%s", junk );      
    fscanf( fp, "%f", &fjunk ); 
	fprintf( fp1, "  %f",fjunk);

    fscanf( fp, "%f", &fjunk ); 
	fprintf( fp1, "  %f",fjunk);

	fgets( junk, 40, fp ); 
	fgets( junk, 40, fp ); 

    fscanf( fp, "%f", &fjunk ); 
	fprintf( fp1, "  %f",fjunk);

    fscanf( fp, "%f", &fjunk ); 
	fprintf( fp1, "  %f",fjunk);

	fprintf( fp1,"\n");
	fclose( fp );
	fclose( fp1 );
}

void checkexpirationdate(void)
{
char cname[80],cname1[40];
int istat;
int year,days;
FILE * pFile;

    sprintf( cname1,"%s.checkthedate",namDynMolBase);


/********************Check the year****************/   

  sprintf( cname,"date +%%Y > %s",cname1);

  istat =  system(cname); 
  if(istat != 0) goto Expire;      
  if ((pFile = fopen (cname1, "r")) == NULL){
                goto Expire;
	} 
  
  fscanf (pFile,"%d",&year);   
//  printf("year is %d\n",year);
  fclose(pFile);
  sprintf( cname,"rm -f %s",cname1);
  istat =  system(cname); 
  if (year > 2009) goto Expire;
  
/*************************************************/   

/************************Check the number of days********************/   
   sprintf( cname,"date +%%j > %s",cname1); 
//   printf("%s\n",cname);
   istat =  system(cname);      
   if(istat != 0) goto Expire; 
   if ((pFile = fopen (cname1, "r")) == NULL){
                goto Expire;
   } 
   fscanf (pFile,"%d",&days);   
   printf("This evaluation version will stop working in %d days\n",(366-days));
   sprintf( cname,"rm -f %s",cname1);
   istat =  system(cname); 
/*********************************************************************/   
  return;   
  Expire:;
  printf("Your evaluation license is expired. Exiting...\n");
  //Exiting();
  exit( 0 );
} 
