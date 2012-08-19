#include "TrigFTKSim/FTKSetup.h"
#include "TrigFTKSim/RoadFinder.h"
#include "TrigFTKSim/FTK_AMBank.h"
#include "TrigFTKSim/tsp/FTKTSPBank.h"
#include "TrigFTKSim/FTK_RawInput.h"
#include "TrigFTKSim/FTKRoadFileOutput.h"
#include "TrigFTKSim/FTKSectorMap.h"
#include "TrigFTKSim/FTKConstantBank.h"
#include "common_fcn.h"
#include "TSystem.h"

#include "boost/lexical_cast.hpp"

#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl_monte_ftk.h>
#include <gsl/gsl_monte_vegas.h>
#include <gsl/gsl_blas.h>

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <ctime>


// Parameters for pattern selection

#define PATT_BASE 0          // The number of patterns to integrate
#define PATT_OFFSET 500          // The number of patterns to integrate

// Just for the example print 
#define EXAMPLE_PATTERN 0 
#define EXAMPLE_SECTOR 48                // The sector selected for example printout in the beginning 


// If activated it uses gsl for matrix-vector multiplication  
//#define USE_GSL 


// Monte Carlo Related patameters
#define MC_METHOD PLAIN		 // 3 available : {FTK, PLAIN, MISER, VEGAS}

#define MIN_CALLS 400000
#define MAX_CALLS 10000000
#define MIN_ERR 0.40


//
//
//     DEBUG PRINTING MODES
//
//
//

// Prints the value of chi2 and exp
//#define VERBOSE_PATTERN


// It activates a printing mode that performs both plain and ftk integration, and prints data to compare them
//#define COMPARE_METHODS 

// Prints every sample and the function value
//#define STUDY_OUTPUT		 	






// @TODO togliere
using namespace std;

// FTKSIM Declarations
// instance of the RoadFinder object
RoadFinder rfobj;

// ss offset (set through config file)
double ss_offset_fraction = 0;

// set the use of the TSP DB banks
int useTSPBank(0);
// if the TSP bank is used it is possible to ask a specific level
int BankLevel(1);
// set if the TSP simulation is performed, 0 not used
int doTSPSim(0);
// minimum TSP coverage, if >1 the way the AM bank is built change
int minTSPCoverage(1);
// store the variable if the TSP bank cache has to be saved
int doMakeCache(0);
// path of the output file
string CachePath("bankcache.root");

int SaveAllRoads(0);

// Plane map pointer, to be set during the initialization
FTKPlaneMap *pmap(0);
FTKPlaneMap *pmap_unused(0);
// Super-strip map, this object transform an hit in a SS index
FTKSSMap *ssmap(0);
FTKSSMap *ssmap_unused(0);
FTKSSMap *ssmap_tsp(0);
// require presence of first layer (pixel B-layer)
bool require_first(false);
// use hashmap internally in AM simulation (allows super-small SS widths)
bool force_am_hashmap(false);
// region map used in this session
FTKRegionMap *rmap(0);
FTKRegionMap *rmap_unused(0);
int CUR_REGION(-1);
// lists containing output of 8L run
std::string scttrk_tracklist;
std::string scttrk_roadlist;

std::string m_badmap_path; 
std::string m_badmap_path2; 

// 4L,8L -> 11L lookup map
FTKSectorMap *scttrk_sectormap(0);

int 	read_commands();
double 	expChi2f (double *x, size_t dim, void *params);
double 	chi2f (double *k, size_t dim, void *params);
double 	gsl_expChi2f (double *x, size_t dim, void *params);
double 	gsl_chi2f (double *k, size_t dim, void *params);
bool 	extract_intervals( const FTKPattern * patt, double *low_vector, double *high_vector, bool *isY  );  
bool 	print_intervals( double *low_vector, double *high_vector, bool *isY  ) ;
bool 	print_intervals_short( double *low_vector, double *high_vector, bool *isY  ) ;
bool 	print_details( FTKPattern * patt, FTK_AMBank * pattBank );
void 	print_k( float ** local_kernel );

// parameters for chi2f
int nconstr, ncoords;

// kernel
gsl_matrix_view kernel_matrix ; // for gsl version
gsl_vector_view kaverages ; // for gsl version
float *kaverage;  
float **kernel;  


/** main function */
int main(int argc, char *argv[]) {
  
  cout << fixed << std::setprecision (3) << "Pattern Volume Calculation : Start." << endl;
 
  // Reading the input
  // 
  FTKSetup &ftkset = FTKSetup::getFTKSetup();
  ftkset.setVerbosity(0);
   
  common_fcn_init(argc, argv);

  // Preparing the output module
  FTKRoadFileOutput *road_output = new FTKRoadFileOutput();
  rfobj.setRoadOutputModule(road_output);
  

  // preliminary parsing of the input file
  ifstream optionfile;
  streambuf *oldrdbuf = cin.rdbuf();


  if (argc > 1) {

    // input option in an external file
    optionfile.open(argv[1]);

    if (!optionfile) {
      cerr << "*** Error reading input file: " << argv[1] << endl;
      return -1;
    }
 
    // input option parsed from stdin
    cin.rdbuf(optionfile.rdbuf());
  }


  // call the parser
  if (read_commands() < 0) {

    cerr << "*** Error parsing the options" << endl;
    return -1;

  }

  cin.rdbuf(oldrdbuf);

  cout << "Init ConstantBank" << endl;
  //ssmap and pmap are global and are initialized in read_commands()
  FTKConstantBank *cnstBank = new FTKConstantBank(pmap->getTotalDim(), "input/data/gcon/corrgen_raw_7L_16M_reg0_sub0.gcon.bz2");
  
  cout << "Get PattBank" << endl;
  //I read the bank id from the conf file but at some point we are going to have read_commands save this somewhere
  FTK_AMBank *pattBank = rfobj.getAMBank(0);

  // Printing the kernel matrix and averages of an example pattern
  FTKPattern myPatt = pattBank->getPattern(EXAMPLE_PATTERN);
   
  print_details( &myPatt, pattBank );
  //
  //  return 0;

  //setting up chi2f parameters 
  nconstr  = cnstBank->getNConstr();
  ncoords  = cnstBank->getNCoords();
  kaverage = cnstBank->getKaverage(EXAMPLE_SECTOR);
  kernel   = cnstBank->getKernel(EXAMPLE_SECTOR);

  ////////// Printing the Kernel matrix and averages
  //
  cout << "\n\n********** kaverage[" << nconstr << "]: ";

  for(int i=0;i<nconstr;++i){
  
	cout << std::setprecision (2)  << kaverage[i] << "\t";

  }

  std::cout << std::endl;

  cout << "\n\n********** kernel["<< nconstr << "][" << ncoords << "] for sector " << EXAMPLE_SECTOR << ": " << endl;

  for(int i=0;i<nconstr;++i){

	for(int j=0;j<ncoords;++j) {
      	
		std::cout << std::setprecision (2) << kernel[i][j] << "\t";

    	}
    
	std::cout << std::endl;
  }


  //std::ofstream pattOutFile;
  //pattOutFile.open ("patt_prob.dat");

  // Sectorized is a Map : int -> list <int>
  // - every key is a sector
  // - every value is the list of pattern indexes that belong to that sector  
  std::map<int, std::list<int> > sectorized;

  std::cout << "\n########## Separate Sectors \n" << std::endl ;
  
  unsigned int limit = PATT_BASE + PATT_OFFSET;

  for ( unsigned int i = PATT_BASE; i < limit ; i ++ ){

	  FTKPattern iterationPatt = pattBank->getPattern(i);      

	  // I add he index of the pattern into the list correspondent to the sector
	  sectorized[boost::lexical_cast<int>(iterationPatt.getSectorID())].push_front(i); 	  

	  if ( i % 1000000 == 0 )
	  {
		  std::cerr << "Separated " << i << " patterns." << std::endl;
	  }

#ifdef VERBOSE_PATTERN
	  std::cout << "\n" <<  i << std::endl;
#endif 
  }


  std::cout << "\n########## Integrating \n" << std::endl ;


  /////////// The integration operation is repeated for every pattern
  // 
  unsigned int 	counter = 0;        /**< counter of the patterns for which the probablility has been calculated */
  unsigned int 	same_counter = 0;	/**< count. of patterns encountered that are from the same sector(usually not used) */
  double 	one_every = 3;			/**< Used if when we want to sample one every N */
  unsigned long int index = 0;		/**< Index of the iteration */
  size_t	calls = 500000;			/**< Number of calls for the monte carlo */ 
  enum		MC_METHOD_TYPE { PLAIN, FTK, MISER, VEGAS}; /**< Monte Carlo Methods */
  int		method = MC_METHOD;		/**< Method that will be used */ 
  unsigned int 	selectedPatterns[PATT_OFFSET];


  // Declaring Monte Carlo Related variables 
  const gsl_rng_type *T;  	// Randomizer type
  gsl_rng *r;			// Randomizer
  double res, err;		// Result and Error 

  // Setting up randomization
  gsl_rng_env_setup ();
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);


  gsl_monte_plain_state *s_plain = gsl_monte_plain_alloc (10);
  gsl_monte_ftk_state 	*s_ftk ; //= gsl_monte_ftk_alloc (10);
//  gsl_monte_ftk_state 	*s_ftk = gsl_monte_ftk_alloc (10);
  gsl_monte_miser_state *s_miser = gsl_monte_miser_alloc (10);
  gsl_monte_vegas_state *s_vegas = gsl_monte_vegas_alloc (10);

  // Performance tuning @TODO togliere
  time_t tstart, tend;

  // 	
  double low[10], high[10];
  bool isY[10];

  tstart = time(0);

#ifdef VERBOSE_PATTERN
  std::cout << "\nStarting the iteration with calls =  " <<  calls << std::endl << std::endl;
#endif 

  std::cout << "index\tcov\tprob\terr\n";

  std::map<int, std::list<int> >::iterator iterator;

  for( iterator = sectorized.begin(); iterator != sectorized.end(); iterator++) 
  {

#ifdef USE_GSL	     
		  double * local_kernel   = (double *)calloc (50, sizeof(double));  
		  float ** tempKernel = cnstBank->getKernel(iterator->first);

		  for(int i = 0 ; i < 5 ; i++)
				  for (int j=0; j < 10 ; j++)
						  local_kernel[i*10+j] = tempKernel[i][j];


		  double * avg = (double *) calloc(5, sizeof(double));

		  for(int i = 0 ; i < 5 ; i++)
				  avg[i] = cnstBank->getKaverage(iterator->first)[i];


		  kernel_matrix = gsl_matrix_view_array(local_kernel, 5, 10);
		  kaverages 	= gsl_vector_view_array(avg, 5 ); 	

		  //for(int i = 0 ; i < 5 ; i++)
		  //	std::cout << "3 " << gsl_vector_get( &kaverages.vector, i) <<  std::endl;

		  gsl_monte_function G = { &gsl_expChi2f, 10, local_kernel };

#else 

		  kaverage = cnstBank->getKaverage(iterator->first);
		  kernel   = cnstBank->getKernel(iterator->first);

		  gsl_monte_function G = { &expChi2f, 10, 0 };
#endif 

	  // Iterating on the patterns in the sector
	  while( ! iterator->second.empty() )
	  {	

		  std::list<int> & myList = iterator->second;

		  int index = myList.front();
		  myList.pop_front();

		  FTKPattern iterationPatt = pattBank->getPattern(index);      

#ifdef VERBOSE_PATTERN
		  std::cout << "\n" << " - Pattern " << index << " : " << iterationPatt << "\n" << std::endl;
#endif 

		  //selectedPatterns[counter] = index; // Saving the selected patterns indexes  



		  //////////// Extracting and printing the intarvals
		  extract_intervals( &iterationPatt, low, high, isY );  

#ifdef STUDY_OUTPUT
		  std::cout << fixed <<"\nNew Pattern : " << index << "\t" << iterationPatt ; 
		  print_intervals_short( low, high, isY  ) ;
#endif 
		  
#ifdef VERBOSE_PATTERN
		  print_intervals( low, high, isY );
#endif
		  //
		  ////////////////////



		  //////// Monte Carlo Integration 
		  //
		  int tempCalls = calls;
		  size_t required_calls = 0;
		  size_t non_zeros = 0;
		  err = res;
		  double err_ratio = err/res;


		  if ( method == FTK )
		  {
			  // Allocating, perfoming integration, and freeing

			  gsl_monte_ftk_integrate_2 (&G, low, high, 10, MIN_CALLS, MAX_CALLS, MIN_ERR, r, s_ftk, &required_calls, &non_zeros,  &res, &err);
			  


#ifdef COMPARE_METHODS


			  double temp_res;
			  double temp_err;
			  size_t temp_non_zeros = 0 ; 

			  gsl_monte_ftk_integrate (&G, low, high, 10, MIN_CALLS, MAX_CALLS, MIN_ERR,  r, s_ftk, &temp_non_zeros, &temp_res, &temp_err);

			  printf("Total Non-zeros_ftk : %d ", (int)non_zeros );
			  printf("Total Non-zeros_plain : %d ",(int) temp_non_zeros );
			 
			  if ( temp_res > res )
			  {
			  std::cout << "\n\nres_diff: " << 1 - ( temp_res / res) << " err_diff " 
					  	<< 1 - ( temp_err / err) << std::endl; 

			  }
			  else
			  {
			  std::cout << "\n\nres_diff: " << 1 - ( res / temp_res ) << " err_diff " 
					  	<< 1 - ( err / temp_err ) << std::endl; 
			  
			  }	  



			  std::cout 	<< index << "\t" << iterationPatt << "\t" << std::setprecision (5) 
					  << scientific << temp_res << "\t" << temp_err <<" \t" << fixed <<  std::setprecision (2) 
					  << err/res << std::setprecision(1) << " \t" << 10000000 
					  << std::setprecision(2) << std::endl ;

#endif

		  }
		  else
		  {	 

			  while( err_ratio > MIN_ERR )
			  {
				  //if ( err_ratio != 1 )
				  //	  std::cout << "Repeating : "	<< "\t" << err_ratio << "\n";

				  switch ( method ){

					  case PLAIN:
						  {
							  // Allocating, perfoming integration, and freeing
							  gsl_monte_plain_integrate (&G, low, high, 10, tempCalls, r, s_plain, &res, &err);
							  break;
						  }
					  case MISER:
						  {
							  gsl_monte_miser_integrate (&G, low, high, 10, tempCalls, r, s_miser, &res, &err);
							  break;	 
						  }	
					  case VEGAS:
						  {
							  gsl_monte_vegas_integrate (&G, low, high, 10, 10000, r, s_vegas, &res, &err);// Warm up
							  gsl_monte_vegas_integrate (&G, low, high, 10, tempCalls, r, s_vegas, &res, &err);
							  break;
						  }	     
				  } // end of switch

				  tempCalls *= 3 ;

				  err_ratio = err/res;

			  } // end of ratio while

		  } // end of else (not FTK method)

		  std::cout 	<< index << "\t" << iterationPatt << "\t" << std::setprecision (5) 
			  	<< scientific << res << "\t" << err <<" \t" << fixed <<  std::setprecision (2) 
			  	<< err/res << std::setprecision(1) << " \t" << required_calls 
			  	<< std::setprecision(2) << std::endl ;

#ifdef VERBOSE_PATTERN
		  // Printing results 
		  std::cout 	<< "At low interval:\n"
			  << std::setprecision (15) << "chi^2 : " << fixed << chi2f( low, 10, 0 ) 
			  << "\nexp(-chi^2/2): " << expChi2f(low, 10, 0) << std::endl  << std::endl ;	

		  std::cout	<< "Integration \n" << std::setprecision (5) << "Value : " 
			  << scientific << res << "\nErr : " << err << fixed << "\n\n";
#endif
		  //
		  ////////////// End of Monte carlo

		  if ( counter % 1000000 == 0 )
		  {
			  std::cerr << "Calculated " << counter << " patterns." << std::endl;
		  }

		  counter++;

	  } // end of list iteration 

#ifdef USE_GSL
	 if(local_kernel)
	  	free( local_kernel ) ;

	 if(avg)
		free( avg ) ;

#endif

  } // end of map iteration 
  //
  /////////////

  std::cout << " End of selected patterns. " << "\n\n";

  // Freeing the monte carlo
  gsl_monte_plain_free (s_plain);
  //gsl_monte_ftk_free (s_ftk);
  gsl_monte_miser_free (s_miser);
  gsl_monte_vegas_free (s_vegas);

  //pattOutFile.close();

  tend = time(0);

  int diff = difftime(tend, tstart) ; 
  double ratio = (double) diff / counter; 

  std::cout 	<< "Calculated  " << counter << " integrals in " << diff << " seconds. \nAvg. " 
	  << ratio << " seconds/pattern" << std::endl;

  return 0;

}



bool print_details( FTKPattern * patt, FTK_AMBank * pattBank ){

  if( ! patt ) // If patt is null
  {
	return false;
  }

  cout << "Printing details of Pattern " <<  *patt << endl;

  cout << " sector coverage " 
       << pattBank->getSectorCoverage(patt->getSectorID())
       << " sectortotal coverage " 
       << pattBank->getTotalCoverage()
       << endl;

  double low[10], high[10];
  bool isY[10];

  int nplanes = patt->getNPlanes();

  for (int iplane=0;iplane!=nplanes;++iplane) { // loop over the planes

    int SSid, section, phioff, phimod, localX, etaoff, etamod, localY, etass, phiss, ix, iy;

    SSid = patt->getSSID(iplane);

    section = 0;//patt.getSectorID();

    cout << "SSid " << SSid 
	 << " plane " << iplane
	 << " section " << section << endl;

    FTKHit hit;
    hit.setPlane(iplane);
    hit.setSector(section);
    phiss = ssmap->getSSPhiWidth(hit); 
    cout << " phiSS " << phiss << endl;
    
    ix = pmap->getDim(iplane,0);
    // indici del vettore delle cooridinate
    iy = pmap->getDim(iplane,1);
    
    if ( ssmap->getSSEtaWidth(hit) != -1 ) // If is pixel
    {
    	etass = ssmap->getSSEtaWidth(hit); 
    	cout << " etaSS " << etass << endl;
    
    	ssmap->decodeSS(SSid, iplane, section, phioff, phimod, localX, etaoff, etamod, localY ); 
    	cout << " phioff " << phioff 
	 << " phimod " << phimod 
	 << " localX " << localX 
	 << " etaoff " << etaoff 
	 << " etamod " << etamod 
	 << " localY " << localY << endl; 
    
    	low[iy] = localY ;
    	high[iy] = localY + etass ;
    	
        isY[ix] = false;
    	isY[iy] = true;
    }
    else  // If is SCT	
    {
		
	ssmap->decodeSS(SSid, iplane, section, phioff, phimod, localX, etaoff); 
    	cout << " phioff " << phioff 
	 << " phimod " << phimod 
	 << " localX " << localX 
	 << " etaoff " << etaoff << endl;
       
        isY[ix] = false;
    }

    low[ix] = localX;
    high[ix] = localX + phiss;   

  }

  print_intervals( low, high,isY  );


}	 

bool print_intervals_short( double *low_vector, double *high_vector, bool *isY  )  {

  std::cout << "Intervals" << std::endl; 


  for (unsigned int i = 0 ; i < 10; i++)
  {  
      if( ! isY[i] )
      {
      	std::cout << " " << low_vector[i] << "-x" << i   << "-" << high_vector[(i)] << " ";  
      }
      else
      {
        std::cout << " " << low_vector[i] << "-y" << i   << "-" << high_vector[(i)] << " ";  
      }
	
      std::cout << "\n";	

  }

}	 

bool print_intervals( double *low_vector, double *high_vector, bool *isY  )  {

  std::cout << "Intervals" << std::endl; 


  for (unsigned int i = 0 ; i < 10; i++)
  {  
      if( ! isY[i] )
      {
      	std::cout << " " << low_vector[i] << "\t< x" << i   << " <\t" << high_vector[(i)] << std::endl;  
      }
      else
      {
        std::cout << " " << low_vector[i] << "\t< y" << i   << " <\t" << high_vector[(i)] << std::endl;  
      }
  }

}	 

bool extract_intervals( const FTKPattern * patt, double *low_vector, double *high_vector, bool *isY  )  {

  /////////// Iterating the planes 
  //
  int nplanes = patt->getNPlanes();

  for (int iplane=0 ; iplane!=nplanes ; ++iplane) { // loop over the planes

    int SSid, section, phioff, phimod, localX, etaoff, etamod, localY, etass, phiss, ix, iy;

    SSid = patt->getSSID(iplane);

    section = 0;//patt.getSectorID();

    FTKHit hit;
    hit.setPlane(iplane);
    hit.setSector(section);
    phiss = ssmap->getSSPhiWidth(hit); 
    
    ix = pmap->getDim(iplane,0);
    iy = pmap->getDim(iplane,1);
    
    if ( ssmap->getSSEtaWidth(hit) != -1 ) // If is pixel
    {
    	etass = ssmap->getSSEtaWidth(hit); 
    
    	ssmap->decodeSS(SSid, iplane, section, phioff, phimod, localX, etaoff, etamod, localY ); 
    
    	low_vector[iy] = localY ;
    	high_vector[iy] = localY + etass ;
    	
        isY[ix] = false ;
    	isY[iy] = true  ;
    }
    else  // If is SCT	
    {
	ssmap->decodeSS(SSid, iplane, section, phioff, phimod, localX, etaoff); 
        isY[ix] = false;
    }

    low_vector[ix] = localX;
    high_vector[ix] = localX + phiss;   

  }

}


void print_k( float ** local_kernel ){

  cout << "\n\n********** kaverage[" << nconstr << "]: ";

  for(int i=0;i<nconstr;++i){
  
	cout << std::setprecision (20)    << local_kernel[i][10] << " ";

  }

  std::cout << std::endl;

  cout << "\n\n********** kernel["<< nconstr << "][" << ncoords << "] for sector " << EXAMPLE_SECTOR << ": " << endl;

  for(int i=0;i<nconstr;++i){

	for(int j=0;j<ncoords;++j) {
      	
		std::cout << std::setprecision (15) << kernel[i][j] << " ";

    	}
    
	std::cout << std::endl;
  }

}


////////////////////// COMMENTED LINES PRESENT IN THE CODE BEFORE JUNE 2012
/////////

  // 	 // decide if a 1D or 2D hit
  // 	 int ndim = iy==-1 ?  1 : 2;
  // 	 if (position[ip]!=endlist[ip]) {
  // 	   if (ndim==1) {
  // 	     newtrk.setCoord(ix,(*position[ip])[0]);
  // 	   }
  // 	   else {
  // 	     newtrk.setCoord(ix,(*position[ip])[0]);
  // 	     newtrk.setCoord(iy,(*position[ip])[1]);
  // 	   }
  // 	 }
  //        }
  //        newtrk.setTrackID(m_comb_id);
  //        newtrk.setHWRejected(HWbase);
  //        newtrk.setHWTrackID(-1);




  //  int totevt = ftkset.getTotEvents();
  //  int ievt_step = totevt > 500 ? totevt / 100 : 1;
  // input module init
  //  cout << "I/O initialization" << endl;
  // instanciate the output module
  //  rfobj.init();
 //  for (int ievt = 0; ievt < totevt; ++ievt) { // events processing loop
 //    if (ievt % ievt_step == 0) {
 //      cout << "Processing evt # " << ievt << " / " << totevt << endl;
 //      if (ievt % (ievt_step * 10) == 0)
 //	ftkset.usageStat();
 //    }
 //    
 //    int res = rfobj.nextEvent();
 //    if (res < 0) { // error
 //      cerr << "*** error reading event # " << ievt << endl;
 //      break;
 //    }
 //  } // events procesing loop
 //
 //  ftkset.usageStat();
 //
  // destroy the I/O modules
  //  delete road_output;
  //  delete rfobj.getDataInputModule();


////////////////
//////////////////////  end of commented lines



/** this function parse the input commands and prepare a simple structure
    used by the main function */
int read_commands() {

	// ftk environemnt
  FTKSetup &ftkset = FTKSetup::getFTKSetup();
  
  const char prompt[] = "PARSE> ";
  
  string line;
  while (replace_getline(cin, line) != 0) {
    if (line[0] == '#' || line.size() == 0) {
     continue;
    }
    
    cout << prompt << line << endl;
    
    // extraxt the pair (key,value)
    istringstream sline(line);
    string key;
    sline >> key;
    if (key == "IGNORE_SIGSEGV") {
      // disable ROOT signal handler (it's buggy and may cause deadlocks after segfault)
      int ignore_sigsegv;
      sline >> ignore_sigsegv;
      if (ignore_sigsegv)
	gSystem->IgnoreSignal(kSigSegmentationViolation, true);
    } else if (key == "events#") {
      int ival;
      sline >> ival;
      cout << prompt << ival << " events will be processed" << endl;
      ftkset.setTotEvents(ival);
    } else if (key == "VERBOSITY") {
      int ival;
      sline >> ival;
      ftkset.setVerbosity(ival);
    } else if (key == "BARREL_ONLY") {
      int ival;
      sline >> ival;
      cout << prompt << "BARREL_ONLY = " << ival << endl;
      ftkset.setBarrelOnly(ival);
    } else if (key == "ENABLE_FTKSIM") {
      int ival;
      sline >> ival;
      cout << prompt << key << " = " << ival << endl;
      ftkset.setEnableFTKSim(ival);
    } else if (key == "SS_FILE") {
      string sval;
      sline >> sval;
      cout << prompt << key << " = " << sval << endl;
      ssmap = new FTKSSMap(rmap, sval.c_str(), force_am_hashmap == false);
    } else if (key == "SS_FILE_UNUSED") {
      string sval;
      sline >> sval;
      cout << prompt << key << " = " << sval << endl;
      rmap_unused = new FTKRegionMap(pmap_unused, rmap->getPath());
      // the map on unused layer skip the boundary check by default
      ssmap_unused = new FTKSSMap(rmap_unused, sval.c_str(), false);
    } else if (key == "SS_FILE_TSP") {
      string sval;
      sline >> sval;
      cout << prompt << key << " = " << sval << endl;
      // the map on unused layer skip the boundary check by default
      ssmap_tsp = new FTKSSMap(rmap, sval.c_str(), false);
    } else if (key == "SS_OFFSET") {
      double sval;
      sline >> sval;
      cout << prompt << key << " = " << sval << endl;
      ss_offset_fraction = sval;
    } else if (key == "MAX_MISSING_PLANES") {
      int ival;
      sline >> ival;
      cout << prompt << key << " = " << ival << endl;
      ftkset.setMaxMissingPlanes(ival);
    } else if (key == "ROADWARRIOR") {
      int ival;
      sline >> ival;
      cout << prompt << key << " = " << ival << endl;
      ftkset.setRoadWarrior(ival);
    } else if (key == "KEEP_REMOVED") {
      int ival;
      sline >> ival;
      ftkset.setKeepRemoved(ival);
    } else if (key == "INPUTFILE") {
      // add a single file to the list
      string sval;
      sline >> sval;
      FTKDataInput *dinput = rfobj.getDataInputModule();
      if (dinput->addFile(sval.c_str()) != -1) {
	cout << "File: " << sval << " added" << endl;
      }
    } else if (key == "FTKDAT_LIST") {
      // add this list of files as input
      string sval;
      sline >> sval;
      FTKDataInput *dinput = rfobj.getDataInputModule();
      int res = dinput->addFilesList(sval.c_str());
      cout << "Added: " << res << " files" << endl;
    } else if (key == "SCTTRK_MODE") {
      int ival;
      sline >> ival;
      ftkset.setSCTtrkMode(ival);
    } else if (key == "SCTTRK_TRACKLIST") {
      sline >> scttrk_tracklist;
    } else if (key == "SCTTRK_ROADLIST") {
      sline >> scttrk_roadlist;
    } else if (key == "SCTTRK_SECTORMAP") {
      string sval;
      sline >> sval;
      // sector map is only needed in SCTtrk mode
      if (ftkset.getSCTtrkMode()) {
	scttrk_sectormap = new FTKSectorMap(sval.c_str());
      }
    } else if (key == "SCTTRK_NLAYERS") {
      int ival;
      sline >> ival;
      rfobj.setSCTtrkNlayers(ival);
    } else if (key == "SCTTRK_NSUBS") {
      int ival;
      sline >> ival;
      rfobj.setSCTtrkNsubs(ival);
    } else if (key == "NSUBREGIONS") {
      int ival;
      sline >> ival;
      rfobj.setNsubregions(ival);
    } else if (key == "SCTTRK_LASTLAYER") {
      int ival;
      sline >> ival;
      rfobj.setSCTtrkLastLayer(ival);
    } else if (key == "SCTTRK_MAX_MISSING_PLANES") {
      int ival;
      sline >> ival;
      rfobj.setSCTtrkMaxMissingPlanes(ival);
    } else if (key == "CUR_REGION") {
      int ival;
      sline >> ival;
      CUR_REGION = ival;
    } else if (key == "CUR_SUBREGION") {
      int ival;
      sline >> ival;
      rfobj.setSubregion(ival);
    } else if (key == "REQUIRE_FIRST") {
      bool ival;
      sline >> ival;
      require_first = ival;
    } else if (key == "FORCE_AM_HASHMAP") {
      bool ival;
      sline >> ival;
      force_am_hashmap = ival;
    } else if (key == "ENCODE_SUBREGION") {
      int ival;
      sline >> ival;
      rfobj.setEncodeSubregion(ival);
    } else if (key == "MAKECACHE") {
      sline >> doMakeCache;
    } else if (key == "CACHEPATH") {
      sline >> CachePath;
    } else if (key == "NBANKS") {
      int ival;
      sline >> ival;
      // set the number of banks
      rfobj.setNBanks(ival);
    } else if (key == "USETSP_BANK") {
      // 0 - not use the TSP ready bank, >0 use the TSP bank
      sline >> useTSPBank;
    } else if (key == "USETSP_DBLEVEL") {
	// choose the level: 
	// 1 - best TSP level,, 2 - 1st group, 3 - 2nd group....
      sline >> BankLevel;
    } else if (key == "USETSP_SIM") {
      /* 0 -  use the TSP bank but only as AM, 
	 1 - build the AM related without further steps
	 2 - build the AM related bank and the DC mask
	 3 - simulate the TSP
       */
      sline >> doTSPSim;
    } else if (key == "USETSP_MINCOVERAGE") {
      /* If this value is >1 the AM is built using all the TSP
	 patterns with the required coverage, than the AM bank is
	 cutted
       */      
      sline >> minTSPCoverage;
    } else if (key == "SAVEALLROADS") {
      sline >> SaveAllRoads;      
    } else if (key == "BADMOD_FILE") {  // set the bad module map path 
      sline >> m_badmap_path;
    } else if (key == "BADMOD_FILE_FOR_HIT") {  // set the bad module map path 
      sline >> m_badmap_path2;
    } else if (key == "BANK") {
      int ibank, maxpatt;
      string bankpath;
      /* retrieve the number of the bank the tha max number 
	 of patterns to read (-1 means all patterns) */
      sline >> ibank >> maxpatt;
      // the next line is the bank path
      replace_getline(cin, bankpath);
      FTK_AMBank *bank = 0x0;
      if (useTSPBank) {
	FTKTSPBank *tspbank = new FTKTSPBank(ibank);
	// set the TSP variables
	tspbank->setSimulateTSP(doTSPSim);
	tspbank->setSSMapTSP(ssmap_tsp);
	tspbank->setTSPMinCoverage(minTSPCoverage);
	tspbank->setMakeCache(doMakeCache>0 ? true : false);
	tspbank->setCachePath(CachePath.c_str());
	bank = tspbank;
      }
      else bank = new FTK_AMBank(ibank);
    
      /* if this is set in a TSP simulation also the roads rejected 
	 by the DC mechanism or by the TSP are recorded. In this case
	 the number of hits is correctly updated to the number of verified
	 layers */
      bank->setSaveAllRoads(SaveAllRoads);
      bank->setRequireFirst(require_first);
      // additional adjustments for SCTtrk mode
      if (ftkset.getSCTtrkMode()) {
	// use STL hashmap version of ss->[patterns] lookup table
	bank->enableLookupMap();
	// always require SCTtrk layer in AM matching
	bank->setRequireLast(true);
      }
      if (force_am_hashmap) {
	bank->enableLookupMap();
      }
      ssmap->setOffsetFraction(ss_offset_fraction);
      bank->setSSMap(ssmap);
      bank->setSSMapUnused(ssmap_unused);

      // set the bad module map path   
      bank->setBadSSMapPath ( m_badmap_path ); 
      bank->setBadSSMapPath2 ( m_badmap_path2 ); 
      
      
      if (bankpath.substr(bankpath.size() - 2).compare("db") == 0) {
	if (bank->readDBBank(bankpath.c_str(), maxpatt, BankLevel) < 0) {
	  delete bank;
	  return -1;
	}
      } else if (bankpath.substr(bankpath.size() - 4) == string("root")) {
	// use the ROOT routine
	if (bank->readROOTBank(bankpath.c_str(), maxpatt) < 0) {
	  // error reading the file
	  delete bank;
	  return -1;
	}
      } else {
	// if the previous check failed try with ASCII file
	if (bank->readASCIIBank(bankpath.c_str(), maxpatt) < 0) {
	  // error reading the file
	  delete bank;
	  return -1;
	}
      }
      rfobj.setAMBank(ibank, bank);
      ftkset.usageStat();
    } else if (key == "RAW_MODE") {
      int ival;
      sline >> ival;
      ftkset.setRawMode(ival);
      FTK_RawInput *imod(0);
      switch (ival) {
      case 1:
	if (!pmap) {
	  cout << "*** RAW data input module needs a valid PMap"
	       << endl;
	  cout << "*** Call PMAP before RAW_MODE" << endl;
	  return -1;
	}
	imod = new FTK_RawInput(pmap);
	rfobj.setDataInputModule(imod);
	break;
      case 2:
	if (!pmap) {
	  cout << "*** RAW data input module needs a valid PMap"
	       << endl;
	  cout << "*** Call PMAP before RAW_MODE" << endl;
	  return -1;
	}
	//pmap_unused = new FTKPlaneMap(pmap->getPath(),true);
	imod = new FTK_RawInput(pmap, pmap_unused);
	imod->setSaveUnused(true);
	rfobj.setDataInputModule(imod);
	rfobj.setSSSearchUnused(true);
	break;
      default:
	cout << "*** RAW_MODE=" << ival << " is not supported" << endl;
	return -1;
      }
    } else if (key == "PMAP_FILE") {
      // load the plane map configuration file
      string sval;
      sline >> sval;
      pmap = new FTKPlaneMap(sval.c_str());
      if (!(*pmap)) {
	cerr << "*** error loading plane map in: " << sval << endl;
	return -1;
      }
    } else if (key == "PMAP_FILE_UNUSED") {
      // load the plane map configuration file
      string sval;
      sline >> sval;
      pmap_unused = new FTKPlaneMap(sval.c_str());
      if (!(*pmap_unused)) {
	cerr << "*** error loading plane map in: " << sval << endl;
	return -1;
      }
    } else if (key == "SCT_PAIR") {
      int ival;
      sline >> ival;
      if (pmap) {
	pmap->setEnablePlanePair(ival);
      } else {
	cerr
	  << "*** impossible to define SCT_PAIR variable, PMAP not instantiated"
	  << endl;
	return -1;
      }
    }
    // First key is for backwards compatibility
    else if (key == "ENABLESCTMISSPAIR" || key == "MAX_MISSING_SCT_PAIRS") {
      int ival;
      sline >> ival;
      ftkset.setMaxMissingSctPairs(ival);
    } else if (key == "RESTRICT_SCT_PAIR_LAYER") {
      int ival;
      sline >> ival;
      ftkset.setRestrictSctPairLayer(ival);
    } else if (key == "RESTRICT_SCT_PAIR_MODULE") {
      int ival;
      sline >> ival;
      ftkset.setRestrictSctPairModule(ival);
      if (ival)
	// This implies RESTRICT_SCT_PAIR_LAYER
	ftkset.setRestrictSctPairLayer(ival);
    } else if (key == "RMAP_FILE") {
      // set rmap file
      string sval;
      sline >> sval;
      rmap = new FTKRegionMap(pmap, sval.c_str());
      if (!(*rmap)) {
	cerr << "*** error using map from: " << sval << endl;
	return -1;
      }
    } else if (key == "DIAG_CLUSTERING") {
      int ival;
      sline >> ival;
      FTKDataInput *dinput = rfobj.getDataInputModule();
      
      // use dynamic casting to guess the input type
      // at the moment is over-shooting
      if (FTK_RawInput *cinput = dynamic_cast<FTK_RawInput*>(dinput)) {
	if (ival > 0) {
	  cinput->setDiagClustering(true);
	} else {
	  cinput->setDiagClustering(false);
	}
      }
    } else if (key == "DUPLICATE_GANGED") {
      int ival;
      sline >> ival;
      FTKDataInput *dinput = rfobj.getDataInputModule();
      
      // use dynamic casting to guess the input type
      // at the moment is over-shooting
      if (FTK_RawInput *cinput = dynamic_cast<FTK_RawInput*>(dinput)) {
	if (ival > 0) {
	  cinput->setDuplicateGanged(true);
	} else {
	  cinput->setDuplicateGanged(false);
	}
      }
    } else if (key == "GANGED_PATTERN_RECOGNITION") {
      int ival;
      sline >> ival;
      FTKDataInput *dinput = rfobj.getDataInputModule();
      
      // use dynamic casting to guess the input type
      // at the moment is over-shooting
      if (FTK_RawInput *cinput = dynamic_cast<FTK_RawInput*>(dinput)) {
	if (ival > 0) {
	  cinput->setGangedPatternRecognition(true);
	} else {
	  cinput->setGangedPatternRecognition(false);
	}
      }
    } else if (key == "SPLIT_BLAYER_MODULES") {
      int ival;
      sline >> ival;
      FTKDataInput *dinput = rfobj.getDataInputModule();
      
      // use dynamic casting to guess the input type
      // at the moment is over-shooting
      if (FTK_RawInput *cinput = dynamic_cast<FTK_RawInput*>(dinput)) {
	if (ival > 0) {
	  cinput->setSplitBlayerModules(true);
	} else {
	  cinput->setSplitBlayerModules(false);
	}
      }
    } else if (key == "MULT_OUT") {
      int ival;
      sline >> ival;
      // set if the output module saves an outfile
      // for each input file, works only with ascii input data
      FTKRoadOutput *outmode = rfobj.getRoadOutputModule();
      if (ival)
	outmode->setMultiOut(true);
      else
	outmode->setMultiOut(false);
    } else if (key == "OUT_FILE") {
      string sval;
      sline >> sval;
      FTKRoadOutput *outmode = rfobj.getRoadOutputModule();
      outmode->setFileName(sval.c_str());
    } else if (key == "OUT_DIR") {
      string sval;
      sline >> sval;
      FTKRoadOutput *outmode = rfobj.getRoadOutputModule();
      outmode->setOutDir(sval.c_str());
    } else {
      cout << prompt << "\"" << key << "\" is not a valid key" << endl;
      continue;
    }
  }
  
  // Initialize 8L SCTtrk road and track input
  if (ftkset.getSCTtrkMode()) // cy
    {
      if (!scttrk_sectormap) {
	cout << "SCTTRK_SECTORMAP was not set in the config file" << endl;
	exit(0);
      }
      rfobj.setSectorMap(scttrk_sectormap);
      FTKDataInput *dinput = rfobj.getDataInputModule();
      if (!rmap) {
	cout << "*** SCTTRK_LIST data input module needs a valid RMap"
	     << endl;
	return -1;
      }
      if (CUR_REGION == -1) {
	cout << "*** Set CUR_REGION in order to load 8L SCTtrk roads"
	     << endl;
	return -1;
      }
      // Sanity check
      if (ftkset.getMaxMissingSctPairs() > rfobj.getSCTtrkMaxMissingPlanes()) {
	cout
	  << "*** WARNING: MAX_MISSING_SCT_PAIRS > SCTTRK_MAX_MISSING_PLANES"
	  << endl << "*** Setting SCTTRK_MAX_MISSING_PLANES = "
	  << ftkset.getMaxMissingSctPairs() << endl;
	rfobj.setSCTtrkMaxMissingPlanes(ftkset.getMaxMissingSctPairs());
      }
      dinput->useRoadsTracks(rmap->getNumRegions(), CUR_REGION);
      int res = dinput->addTrackFilesList(CUR_REGION,
					  scttrk_tracklist.c_str());
      cout << "Added: " << res << " track files" << endl;
      res = dinput->addRoadFilesList(CUR_REGION, scttrk_roadlist.c_str());
      cout << "Added: " << res << " road files" << endl;
      dinput->initRoadsTracks();
    }
  // Sanity check for the default case
  else if (ftkset.getMaxMissingSctPairs() > ftkset.getMaxMissingPlanes()) {
    cout << "*** WARNING: MAX_MISSING_SCT_PAIRS > MAX_MISSING_PLANES"
	 << endl << "*** Setting MAX_MISSING_PLANES = "
	 << ftkset.getMaxMissingSctPairs() << endl;
    ftkset.setMaxMissingPlanes(ftkset.getMaxMissingSctPairs());
  }
  
  return 0;

}





double gsl_chi2f (double *x, size_t dim, void *params){

	double chi2 = 0;

	// preparing result vector
	double results [ncoords];
	gsl_vector_view results_vector = gsl_vector_view_array(results, 5);
	gsl_vector_set_zero (&results_vector.vector);
	
	// Getting input vector
	gsl_vector_view X = gsl_vector_view_array(x, 10);

	// Multiplication. See gsl regerence
	gsl_blas_dgemv ( CblasNoTrans, 1, & kernel_matrix.matrix, & X.vector , 1, &results_vector.vector );
	// adding the averages
	gsl_vector_add ( &results_vector.vector, &kaverages.vector );

	// Calculating the square of the vector
	gsl_vector_mul ( &results_vector.vector , &results_vector.vector );

	// calculating the sum 
	for( int i = 0 ;i<nconstr;i++)
	{
		chi2 += gsl_vector_get ( &results_vector.vector, i);
	}

	return chi2;

}

double gsl_expChi2f (double *x, size_t dim, void *params){

	return exp(-gsl_chi2f(x, dim, params)/2);

}

double chi2f (double *x, size_t dim, void *params){

	double chi2(0);

	for(int i=0;i<nconstr;i++)
	{
		double s = kaverage[i];
		for(int j=0;j<ncoords;j++) {
			s += kernel[i][j]*x[j];
		}
		chi2 += s*s;
	}

	return chi2;

}


double expChi2f (double *x, size_t dim, void *params){

	// This has to be moved inside our function	
	for( unsigned int i = 0 ; i < dim  ; i++  )	
		x[i] = (int) x[i] ; 

	#ifdef STUDY_OUTPUT
		for( unsigned int i = 0 ; i < dim  ; i++  )	
			std::cout << fixed << std::setprecision (0) << x[i] << " "  ; 
		
		double mychi2 = chi2f(x, dim, params);
		

		std::cout << std::setprecision (15) << chi2f(x, dim, params) << " "  ; 
		std::cout << exp(- mychi2/2) << std::setprecision (2) << std::endl  ; 

	#endif 
	return exp(-chi2f(x, dim, params)/2);

}




/*  VERSIONE VERkaverageOSISSIMA
double chi2f (double *x, size_t dim, void *params){

	double chi2(0);

	float ** mykernel = (float **) params;

	printf("%s " , "x : ");
	for(int j=0; j<ncoords; j++) {
		printf("%f " , x[j]);
	}
	printf("\n" );


	print_k(mykernel);

	for(int i=0;i<nconstr;i++)
	{
		double s = mykernel[i][10];
		printf("Inizio : %f\n", mykernel[i][10] );
		for(int j=0;j<ncoords;j++) {
			s += mykernel[i][j]*x[j];
			printf("%f * %f = %f \n", mykernel[i][j] , x[j], s);
		}
		printf("%s %f\n", "s : " , s);
		chi2 += s*s;
	}

	return chi2;
}*/
