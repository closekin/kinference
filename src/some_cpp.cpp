// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(BH)]]

#include <RcppArmadillo.h>
// #include <boost/tuple/tuple.hpp>
#include <queue>

using namespace Rcpp;
using namespace std;


// Assertions: ought to control this via NDEBUG or similar def
#define STRINGIZE(x) STRINGIZE2(x)
#define STRINGIZE2(x) #x
#define LINE_STRING STRINGIZE(__LINE__)
bool stoppity(
    std::string errmsg
){
  Function Rstop = Environment::base_env()[ "stop"];
  Rcpp::StringVector Rerrmsg (1);
  Rerrmsg[0] = errmsg;
  Rstop( Rerrmsg);
  return false;
}
// Rcpp::stop...
#define ASSERTO(EX) (void)((EX) || (stoppity( "Failed: " #EX " in " __FILE__ " line " LINE_STRING "!"), 0))


// [[Rcpp::export]]
SEXP paircomps(
  IntegerMatrix pair_geno, // n_geno, n_geno -> n_genopairs
  NumericMatrix LOD, // n_genopairs, n_loci,
  RawMatrix geno1, // n_loci, n_samps1,
  RawMatrix geno2, // n_loci, n_samps2
  bool symmo,
  int granulum,
  int granulum_loci
  ) {
  BEGIN_RCPP
    // Avoid scoping woes by doing the include-file right here
#if defined( MVBDEBUG)
#include <../genobasics-vistu/r_int_defs_debug.h>
#endif
    int npc; // #comps
    int n_geno;
    int n_samps1;
    int n_samps2;
    int n_loci;
    int n_genopairs;
    int start_samp1;
    int i;
    int j;
    int i_max;
    int j_max;
    int start_loc;
    int iloc;
    int loc_max;
    int ij_out;
    int g1g2;
    double this_PLOD;
    double this_LOD;

    n_geno = pair_geno. nrow();
    ASSERTO( pair_geno.ncol() == n_geno);

    n_genopairs = LOD.nrow(); // allow triangular storage

    // Ensure matching size of LOD and pair_geno
    // Can't seem to directly get length of pair_geno@what...
    // ... even though ALL R objects have a length
    // ... sigh
    Rcpp::CharacterVector whattr( pair_geno. attr( "what"));
    ASSERTO( whattr.size() == n_genopairs);

    n_loci = LOD.ncol();
    ASSERTO( geno1.nrow() == n_loci);
    n_samps1 = geno1.ncol();
    ASSERTO( geno2.nrow() == n_loci);
    n_samps2 = geno2.ncol();

    if( symmo) {
      ASSERTO( n_samps1 == n_samps2);
      npc = (n_samps1 * (n_samps1-1)) / 2;
    } else {
      npc = n_samps1 * n_samps2;
    };

    NumericVector PLOD ( npc); // vector not matrix for efficient triangular symm

    // Attempt at cache efficiency, by blocking--- avoid
    // re-loading within j-loop or loc-loop
    start_loc = 0;
    while(start_loc <= n_loci) {
      loc_max = std::min( n_loci, start_loc + granulum_loci);
      start_samp1 = 0;
      ij_out = 0;
      while( start_samp1 <= n_samps1) {
        i_max = std::min( n_samps1, start_samp1 + granulum);
        for( i = start_samp1; i < i_max; i++){
          j_max = symmo ? i : n_samps2;
          for( j = 0; j < j_max; j++) {
            this_PLOD = 0;
            for( iloc = start_loc; iloc < loc_max; iloc++ ) {
              g1g2 = pair_geno( geno1( iloc, i), geno2( iloc, j));
              ASSERTO( (g1g2 >= 1) && (g1g2 <= n_genopairs));
              this_LOD = LOD( g1g2-1, iloc); // Effing 0-base...
              this_PLOD += this_LOD;
            };

            PLOD( ij_out) += this_LOD;
            ij_out++;
          }; // for j
        }; // for i
        start_samp1 += granulum;
      }; // for istart
      start_loc += granulum_loci;
    }; // for locus block

    return PLOD;
  END_RCPP
};

// [[Rcpp::export]]
SEXP HSP_paircomps_lots(
  IntegerMatrix pair_geno, // n_geno, n_geno -> n_genopairs
  NumericMatrix LOD, // n_genopairs, n_loci,
  RawMatrix geno1, // n_loci, n_samps1,
  RawMatrix geno2, // n_loci, n_samps2
  bool symmo,
  double eta,
  double min_keep_PLOD,
  int keep_n, // number of PLODs to return
  NumericVector bins
  ) {
BEGIN_RCPP
    // Avoid scoping woes by doing the include-file right here
#if defined( MVBDEBUG)
#include <../genobasics-vistu/r_int_defs_debug.h>
#endif

  int n_geno;
  int n_samps1;
  int n_samps2;
  int n_loci;
  int n_genopairs;
  int n_bins;
  int n_kept;
  int n_sub_PLOD; // discarded ones
  int i;
  int j;
  int j_max;
  int iloc;
  int g1g2;
  // int g1; // debugging
  // int g2;  // debugging
  int ibin;
  double this_PLOD;
  double this_LOD;
  double mean_sub_PLOD;
  double mean_sub_PLOD2;
  double var_sub_PLOD;

  n_geno = pair_geno. nrow();
  ASSERTO( pair_geno.ncol() == n_geno);

  n_genopairs = LOD.nrow(); // allow triangular storage

  // Ensure matching size of LOD and pair_geno
  // Can't seem to directly get length of pair_geno@what...
  // ... even though ALL R objects have a length
  // ... sigh
  Rcpp::CharacterVector whattr( pair_geno. attr( "what"));
  ASSERTO( whattr.size() == n_genopairs);

  n_loci = LOD.ncol();
  ASSERTO( geno1.nrow() == n_loci);
  n_samps1 = geno1.ncol();
  ASSERTO( geno2.nrow() == n_loci);
  n_samps2 = geno2.ncol();

  if( symmo) {
    ASSERTO( n_samps1 == n_samps2);
  };

  n_bins = bins.size();
  IntegerVector n_PLODs_below( n_bins);

  // Need to grow the kept-values, so std::vector is apparently better
  // NumericVector big_PLOD( guess_n_keep);
  // IntegerVector big_i( guess_n_keep);
  // IntegerVector big_j( guess_n_keep);
  //
  // Now get rid of these and use containers
  std::vector<double> big_PLOD;
  std::vector<int>  big_i;
  std::vector<int> big_j;
  // Define a new type, which is the PLOD and i and j
  typedef std::tuple<double, int, int> PLODder;
  // need to define a way in which we sort a priority_queue of tuples...
  // https://stackoverflow.com/questions/5712172/how-can-i-store-3-integer-in-priority-queue
  class my_greater  {
  public:
    bool operator() (const PLODder& arg1, const PLODder& arg2) const{
      return get<0>(arg1) > get<0>(arg2);
      return false;
    }
  };
  std::priority_queue< PLODder, std::vector<PLODder>, my_greater> PLODing_along;


  // check compilation; at one point couldn't find scope
  // now can, but following line doesn't work; prolly doesn't matter
  // n_kept = XINTEGER(n_PLODs_below)[1];

  n_kept = 0;
  n_sub_PLOD = 0;
  mean_sub_PLOD = 0;
  mean_sub_PLOD2 = 0;

  for( i = 0; i < n_samps1; i++){
    // Could maybe speed this up a bit...
    // but next doesn't seem to work. Need something that gives a pointer to
    // first-row-of-column-X
    // NumericMatrix::Column geno1_i = geno1(_,i); // google "make your own apply Rcpp"
    // .... pair_geno( geno1_i[ iloc], geno2_j[ iloc]
    // Bigger speedup would come from "blocking" as per mat-mult to reduce
    // cache misses

    // IntegerMatrix::Column geno1_i = geno1(_,i);
    j_max = symmo ? i : n_samps2; // j_max is i if symmo, or n_samps2 otherwise
    for( j = 0; j < j_max; j++) {
      // IntegerMatrix::Column geno2_j = geno2(_,j);
      // calculate the PLOD
      this_PLOD = 0;
      for( iloc = 0; iloc < n_loci; iloc++ ) {
        g1g2 = pair_geno( geno1( iloc, i)-1, geno2( iloc, j)-1); // Effing 0-base...
        // g1g2 = pair_geno( geno1_i[ iloc]-1, geno2_j[ iloc]-1); // Effing 0-base...
        ASSERTO( (g1g2 >= 1) && (g1g2 <= n_genopairs));
        this_LOD = LOD( g1g2-1, iloc); // Effing 0-base...
        this_PLOD += this_LOD;
      }; // done calculating PLOD


      // now, what to do with it?!

      if( this_PLOD < eta) { // common; predictable branch
        mean_sub_PLOD += this_PLOD;
        mean_sub_PLOD2 += this_PLOD*this_PLOD;
        n_sub_PLOD += 1;
      };

      if( this_PLOD > min_keep_PLOD) { // rare; predictable branch
        // build our new friend to insert
        PLODder new_PLOD_on_the_block = PLODder(this_PLOD, i, j);

        // if we didn't fill up yet, then just push away
        if( PLODing_along.size() < keep_n){
          PLODing_along.push(new_PLOD_on_the_block);
        }else{
          // uh oh we filled-up! check to see if we really need to add
          // this one (is it better than what's in there?), then exile
          // the last element to the garbage and induct our new best friend
          if( get<0>(PLODing_along.top()) < this_PLOD){
            PLODing_along.pop();
            PLODing_along.push(new_PLOD_on_the_block);
          }
        }
      };

      // how did that effect this bit vvvvv
      // this is where we should update re: issue #36.
      for( ibin = 0; ibin < n_bins; ibin++) { // avoid if() which is slow
        n_PLODs_below( ibin) += (int) ( this_PLOD < bins( ibin));
      };
    }; // for j
  }; // for i

  mean_sub_PLOD /= n_sub_PLOD;
  var_sub_PLOD = mean_sub_PLOD2 / n_sub_PLOD - mean_sub_PLOD*mean_sub_PLOD;

  // n_PLODs_below is cumul but should really be n_PLODs_in_bin, so...
  for( ibin = n_bins; ibin > 0; ibin--){
    n_PLODs_below[ ibin] -= n_PLODs_below[ ibin-1];
  };

  // build big_i, big_j, big_PLOD
  big_i.reserve(PLODing_along.size());
  big_j.reserve(PLODing_along.size());
  big_PLOD.reserve(PLODing_along.size());

  // unwrap that nice structure I guess :'(
  while (!PLODing_along.empty()) {
    big_PLOD.push_back(get<0>(PLODing_along.top())); // "effing 0-base" is what I should write here :P
    big_i.push_back(get<1>(PLODing_along.top()) + 1);
    big_j.push_back(get<2>(PLODing_along.top()) + 1);
    PLODing_along.pop();
  }


  // make a list object using wrap() for the stdvectors
return( Rcpp::List::create(
      Rcpp::Named("big_PLOD")=wrap( big_PLOD),
      Rcpp::Named("big_i")=wrap( big_i),
      Rcpp::Named("big_j")=wrap( big_j),
      Rcpp::Named( "n_PLODs_in_bin")=n_PLODs_below,
      Rcpp::Named( "mean_sub_PLOD")=mean_sub_PLOD,
      Rcpp::Named( "var_sub_PLOD")=var_sub_PLOD
  ));

END_RCPP
};




// [[Rcpp::export]]
SEXP POP_paircomps_lots(
  RawMatrix geno1, // n_loci, n_samps1,
  RawMatrix geno2, // n_loci, n_samps2
  bool symmo,
  double eta,
  int max_keep_Nexclu,
  int keep_n, // number of pairs to return
  NumericVector bins,
  int AAO,
  int BBO
  ) {
BEGIN_RCPP
    // Avoid scoping woes by doing the include-file right here
#if defined( MVBDEBUG)
#include <../genobasics-vistu/r_int_defs_debug.h>
#endif

  int n_samps1;
  int n_samps2;
  int n_loci;
  int n_bins;
  int n_kept;
  int n_hi_Nexclu; // discarded ones
  int i;
  int j;
  int j_max;
  int g2; // debugging
  int iloc;
  int iiloc;
  int n_AAOs1_i;
  int n_BBOs1_i;
  int ibin;
  int this_Nexclu;
  double mean_hi_Nexclu;
  double mean_hi_Nexclu2;
  double var_hi_Nexclu;

  // Ensure matching size of LOD and pair_geno
  // Can't seem to directly get length of pair_geno@what...
  // ... even though ALL R objects have a length
  // ... sigh

  n_loci = geno1.nrow();
  n_samps1 = geno1.ncol();
  ASSERTO( geno2.nrow() == n_loci);
  n_samps2 = geno2.ncol();

  if( symmo) {
    ASSERTO( n_samps1 == n_samps2);
  };

  IntegerVector which_AAOs1_i( n_loci);
  IntegerVector which_BBOs1_i( n_loci);

  n_bins = bins.size();
  IntegerVector n_Nexclu_below( n_bins);

  // Need to grow the kept-values, so std::vector is apparently better
  //std::vector<double> big_Nexclu;
  //std::vector<int>  big_i;
  //std::vector<int> big_j;
  // Now get rid of these and use containers
  std::vector<double> big_Nexclu;
  std::vector<int>  big_i;
  std::vector<int> big_j;
  // Define a new type, which is the Nexclu and i and j
  typedef std::tuple<double, int, int> Nexcluder;
  // define what "better" means
  class my_less  {
  public:
    bool operator() (const Nexcluder& arg1, const Nexcluder& arg2) const{
      return get<0>(arg1) <= get<0>(arg2);
      return false;
    }
  };
  std::priority_queue< Nexcluder, std::vector<Nexcluder>, my_less> Nexcluing_along;

  // check compilation; at one point couldn't find scope
  // now can, but following line doesn't work; prolly doesn't matter
  // n_kept = XINTEGER(n_PLODs_below)[1];

  n_kept = 0;
  n_hi_Nexclu = 0;
  mean_hi_Nexclu= 0;
  mean_hi_Nexclu2 = 0;

  for( i = 0; i < n_samps1; i++){
    // Could maybe speed this up a bit...
    // but next doesn't seem to work. Need something that gives a pointer to first-row-of-column-X

    // For speed, only check exlcus over loci where g1 == AAO (then just check if g2 == BBO there) and conversely
    // Set up info for i; only done once

    n_AAOs1_i = 0;
    n_BBOs1_i = 0;
    for( iloc = 0; iloc < n_loci; iloc++) {
      if( geno1( iloc, i) == AAO) {
        which_AAOs1_i[ n_AAOs1_i] = iloc;
        n_AAOs1_i += 1;
      } else if( geno1( iloc, i) == BBO) {
        which_BBOs1_i[ n_BBOs1_i] = iloc;
        n_BBOs1_i += 1;
      }
    };

    j_max = symmo ? i : n_samps2;
    for( j = 0; j < j_max; j++) {
      this_Nexclu = 0;
      for( iiloc = 0; iiloc < n_AAOs1_i; iiloc++) {
        iloc = which_AAOs1_i[ iiloc];
        g2 = geno2( iloc, j);
        this_Nexclu += (int)( geno2( which_AAOs1_i[ iiloc], j) == BBO);
      };

      for( iiloc = 0; iiloc < n_BBOs1_i; iiloc++) {
        this_Nexclu += (int)( geno2( which_BBOs1_i[ iiloc], j) == AAO);
      };

      if( this_Nexclu > eta) { // common; predictable branch
        mean_hi_Nexclu += this_Nexclu;
        mean_hi_Nexclu2 += this_Nexclu*this_Nexclu;
        n_hi_Nexclu += 1;
      };

      if( this_Nexclu <= max_keep_Nexclu) { // rare; predictable branch
        //big_Nexclu. push_back( this_Nexclu);
        //big_i. push_back( i+1); // Effing 0-base...
        //big_j. push_back( j+1); // Effing 0-base...

        Nexcluder new_Nexclu_on_the_block = Nexcluder(this_Nexclu, i, j);

        // if we didn't fill up yet, then just push away
        if( Nexcluing_along.size() < keep_n){
          Nexcluing_along.push(new_Nexclu_on_the_block);
        }else{
          // uh oh we filled-up! check to see if we really need to add
          // this one (is it better than what's in there?), then exile
          // the last element to the garbage and induct our new best friend
          if( get<0>(Nexcluing_along.top()) > this_Nexclu){
            Nexcluing_along.pop();
            Nexcluing_along.push(new_Nexclu_on_the_block);
          }
        }

        // as before
        n_kept += 1;
      };

      // update this for loop for issue #36 with the new even-spaced bins.
      for( ibin = 0; ibin < n_bins; ibin++) { // avoid if() which is slow
        n_Nexclu_below( ibin) += (int) ( this_Nexclu < bins( ibin));
      };
    }; // for j
  }; // for i

  mean_hi_Nexclu /= n_hi_Nexclu;
  var_hi_Nexclu = mean_hi_Nexclu2 / n_hi_Nexclu - mean_hi_Nexclu*mean_hi_Nexclu;

  // n_Nexclu_below is cumul but should really be n_Nexclu_in_bin, so...
  for( ibin = n_bins; ibin > 0; ibin--){
    n_Nexclu_below[ ibin] -= n_Nexclu_below[ ibin-1];
  };


  // build big_i, big_j, big_Nexclu
  big_i.reserve(Nexcluing_along.size());
  big_j.reserve(Nexcluing_along.size());
  big_Nexclu.reserve(Nexcluing_along.size());

  // unwrap that nice structure I guess :'(
  while (!Nexcluing_along.empty()) {
    big_Nexclu.push_back(get<0>(Nexcluing_along.top())); // "effing 0-base" is what I should write here :P
    big_i.push_back(get<1>(Nexcluing_along.top()) + 1);
    big_j.push_back(get<2>(Nexcluing_along.top()) + 1);
    Nexcluing_along.pop();
  }



  // make a list object using wrap() for the stdvectors
return( Rcpp::List::create(
      Rcpp::Named("big_Nexclu")=wrap( big_Nexclu),
      Rcpp::Named("big_i")=wrap( big_i),
      Rcpp::Named("big_j")=wrap( big_j),
      Rcpp::Named( "n_Nexclu_in_bin")=n_Nexclu_below,
      Rcpp::Named( "mean_hi_Nexclu")=mean_hi_Nexclu,
      Rcpp::Named( "var_hi_Nexclu")=var_hi_Nexclu
  ));

END_RCPP
};

// [[Rcpp::export]]
SEXP POP_wt_paircomps_lots(
  RawMatrix geno1, // n_loci, n_samps1,
  RawMatrix geno2, // n_loci, n_samps2
  NumericVector w, // n_loci
  bool symmo,
  double eta,
  double max_keep_wpsex,
  int keep_n, // number of pairs to return
  NumericVector bins,
  int AAO,
  int BBO
) {
  BEGIN_RCPP
    // Avoid scoping woes by doing the include-file right here
#if defined( MVBDEBUG)
#include <../genobasics-vistu/r_int_defs_debug.h>
#endif

  int n_samps1;
  int n_samps2;
  int n_loci;
  int n_bins;
  int n_kept;
  int n_hi_wpsex; // discarded ones
  int i;
  int j;
  int j_max;
  int g2; // debugging
  int iloc;
  int iiloc;
  int n_AAOs1_i;
  int n_BBOs1_i;
  int ibin;
  double this_wpsex;
  double mean_hi_wpsex;
  double mean_hi_wpsex2;
  double var_hi_wpsex;

  // Ensure matching size of LOD and pair_geno
  // Can't seem to directly get length of pair_geno@what...
  // ... even though ALL R objects have a length
  // ... sigh

  n_loci = geno1.nrow();
  n_samps1 = geno1.ncol();
  ASSERTO( w. size() == n_loci);
  ASSERTO( geno2.nrow() == n_loci);
  n_samps2 = geno2.ncol();

  if( symmo) {
    ASSERTO( n_samps1 == n_samps2);
  };

  IntegerVector which_AAOs1_i( n_loci);
  IntegerVector which_BBOs1_i( n_loci);

  n_bins = bins.size();
  IntegerVector n_wpsex_below( n_bins);

  // Need to grow the kept-values, so std::vector is apparently better
  //std::vector<double> big_wpsex;
  //std::vector<int>  big_i;
  //std::vector<int> big_j;
  // Now get rid of these and use containers
  std::vector<double> big_WPSEX;
  std::vector<int>  big_i;
  std::vector<int> big_j;
  // Define a new type, which is the WPSEX and i and j
  typedef std::tuple<double, int, int> WPSEXder;
  // define what "better" means
  class my_less  {
  public:
    bool operator() (const WPSEXder& arg1, const WPSEXder& arg2) const{
      return get<0>(arg1) <= get<0>(arg2);
      return false;
    }
  };
  std::priority_queue< WPSEXder, std::vector<WPSEXder>, my_less> WPSEXing_along;

  // check compilation; at one point couldn't find scope
  // now can, but following line doesn't work; prolly doesn't matter
  // n_kept = XINTEGER(n_PLODs_below)[1];

  n_kept = 0;
  n_hi_wpsex = 0;
  mean_hi_wpsex= 0;
  mean_hi_wpsex2 = 0;

  for( i = 0; i < n_samps1; i++){
    // Could maybe speed this up a bit...
    // but next doesn't seem to work. Need something that gives a pointer to first-row-of-column-X

    // For speed, only check exlcus over loci where g1 == AAO (then just check if g2 == BBO there) and conversely
    // Set up info for i; only done once

    n_AAOs1_i = 0;
    n_BBOs1_i = 0;
    for( iloc = 0; iloc < n_loci; iloc++) {
      if( geno1( iloc, i) == AAO) {
        which_AAOs1_i[ n_AAOs1_i] = iloc;
        n_AAOs1_i += 1;
      } else if( geno1( iloc, i) == BBO) {
        which_BBOs1_i[ n_BBOs1_i] = iloc;
        n_BBOs1_i += 1;
      }
    };

    j_max = symmo ? i : n_samps2;
    for( j = 0; j < j_max; j++) {
      this_wpsex = 0;
      for( iiloc = 0; iiloc < n_AAOs1_i; iiloc++) {
        iloc = which_AAOs1_i[ iiloc];
        g2 = geno2( iloc, j);
        this_wpsex += w[ iloc-1] * (int)( geno2( which_AAOs1_i[ iiloc], j) == BBO);
      };

      for( iiloc = 0; iiloc < n_BBOs1_i; iiloc++) {
        iloc = which_AAOs1_i[ iiloc];
        this_wpsex += w[ iloc-1] * (int)( geno2( which_BBOs1_i[ iiloc], j) == AAO);
      };

      if( this_wpsex > eta) { // common; predictable branch
        mean_hi_wpsex += this_wpsex;
        mean_hi_wpsex2 += this_wpsex*this_wpsex;
        n_hi_wpsex += 1;
      };

      if( this_wpsex <= max_keep_wpsex) { // rare; predictable branch
        //big_wpsex. push_back( this_wpsex);
        //big_i. push_back( i+1); // Effing 0-base...
        //big_j. push_back( j+1); // Effing 0-base...
        WPSEXder new_WPSEX_on_the_block = WPSEXder(this_wpsex, i, j);

        // if we didn't fill up yet, then just push away
        if( WPSEXing_along.size() < keep_n){
          WPSEXing_along.push(new_WPSEX_on_the_block);
        }else{
          // uh oh we filled-up! check to see if we really need to add
          // this one (is it better than what's in there?), then exile
          // the last element to the garbage and induct our new best friend
          if( get<0>(WPSEXing_along.top()) > this_wpsex){
            WPSEXing_along.pop();
            WPSEXing_along.push(new_WPSEX_on_the_block);
          }
        }

        // as before
        n_kept += 1;
      };

      for( ibin = 0; ibin < n_bins; ibin++) { // avoid if() which is slow
        n_wpsex_below( ibin) += (int) ( this_wpsex < bins( ibin));
      };
    }; // for j
  }; // for i

  mean_hi_wpsex /= n_hi_wpsex;
  var_hi_wpsex = mean_hi_wpsex2 / n_hi_wpsex - mean_hi_wpsex*mean_hi_wpsex;

  // n_wpsex_below is cumul but should really be n_wpsex_in_bin, so...
  for( ibin = n_bins; ibin > 0; ibin--){
    n_wpsex_below[ ibin] -= n_wpsex_below[ ibin-1];
  };

  // build big_i, big_j, big_Nexclu
  big_i.reserve(WPSEXing_along.size());
  big_j.reserve(WPSEXing_along.size());
  big_WPSEX.reserve(WPSEXing_along.size());

  // unwrap that nice structure I guess :'(
  while (!WPSEXing_along.empty()) {
    big_WPSEX.push_back(get<0>(WPSEXing_along.top())); // I should write "effing 0-base" here :P
    big_i.push_back(get<1>(WPSEXing_along.top()) + 1);
    big_j.push_back(get<2>(WPSEXing_along.top()) + 1);
    WPSEXing_along.pop();
  }

  // make a list object using wrap() for the stdvectors
  return( Rcpp::List::create( 
    Rcpp::Named("big_wpsex")=wrap( big_WPSEX),
    Rcpp::Named("big_i")=wrap( big_i),
    Rcpp::Named("big_j")=wrap( big_j),
    Rcpp::Named( "n_wpsex_in_bin")=n_wpsex_below,
    Rcpp::Named( "mean_hi_wpsex")=mean_hi_wpsex,
    Rcpp::Named( "var_hi_wpsex")=var_hi_wpsex
  ));

  END_RCPP
};


// [[Rcpp::export]]
SEXP DUP_paircomps_lots(
  RawMatrix geno1, // n_loci, n_samps1,
  RawMatrix geno2, // n_loci, n_samps1,
  bool symmo, // should only be true if geno1==geno2
  double max_diff_genos, // max tolerated discrepant 4way genos to be a "true duplicate"
  int keep_n // number of pairs to return
) {
  BEGIN_RCPP
    // Avoid scoping woes by doing the include-file right here
#if defined( MVBDEBUG)
#include <../genobasics-vistu/r_int_defs_debug.h>
#endif

  int n_samps1;
  int n_samps2;
  int n_loci;
  int n_kept;
  int i;
  int j;
  int j_max;
  int iloc;
  int this_ndiff;

  n_loci = geno1.nrow();
  n_samps1 = geno1.ncol();
  ASSERTO( geno2.nrow() == n_loci);
  n_samps2 = geno2.ncol();

  if( symmo) {
    ASSERTO( n_samps1 == n_samps2);
  };


  // Need to grow the kept-values, so std::vector is apparently better
  //std::vector<double> big_similar;
  //std::vector<int>  big_i;
  //std::vector<int> big_j;
  // Now get rid of these and use containers
  std::vector<double> big_similar;
  std::vector<int>  big_i;
  std::vector<int> big_j;
  // Define a new type, which is the similarity and i and j
  typedef std::tuple<double, int, int> similarder;
  // define what "better" means
  class my_less  {
  public:
    bool operator() (const similarder& arg1, const similarder& arg2) const{
      return get<0>(arg1) <= get<0>(arg2);
      return false;
    }
  };
  std::priority_queue< similarder, std::vector<similarder>, my_less> similaring_along;

  // check compilation; at one point couldn't find scope
  // now can, but following line doesn't work; prolly doesn't matter
  // n_kept = XINTEGER(n_PLODs_below)[1];

  n_kept = 0;
  IntegerVector known_dup( n_samps2, 0);

  for( i = 0; i < n_samps1; i++){
    // "Fuzzy hashing" would really speed this up...
    // preprocess each genotype into a much smaller hash that is guaranteed to differ
    // by "a lot" if genos are very different; only check detailed genos
    // when hashes of i and j are "close".

    j_max = symmo ? i : n_samps2;
    for( j = 0; j < j_max; j++) {
      if( !symmo || !known_dup[ j]) {

        this_ndiff = 0;
        for( iloc = 0; iloc < n_loci; iloc++) {
          // Faster version could "unroll" this a bit and only do the if-check every 10 loci or so
          this_ndiff += (int)( geno2( iloc, j) != geno1( iloc, i));
          if( this_ndiff > max_diff_genos) { // not a dup
        break;
          };
        };

        if( this_ndiff > max_diff_genos) {
      continue;
        } else {
          known_dup[ j] = 1;
          //big_similar. push_back( this_ndiff);
          //big_i. push_back( i+1); // Effing 0-base...
          //big_j. push_back( j+1); // Effing 0-base...
          similarder new_similar_on_the_block = similarder(this_ndiff, i, j);

          // if we didn't fill up yet, then just push away
          if( similaring_along.size() < keep_n){
            similaring_along.push(new_similar_on_the_block);
          }else{
            // uh oh we filled-up! check to see if we really need to add
            // this one (is it better than what's in there?), then exile
            // the last element to the garbage and induct our new best friend
            if( get<0>(similaring_along.top()) > this_ndiff){
              similaring_along.pop();
              similaring_along.push(new_similar_on_the_block);
            }
          }

          // as before
          n_kept += 1;
        };
      }; // if j not already known to be a dup
    }; // for j
  }; // for i

  // build big_i, big_j, big_Nexclu
  big_i.reserve(similaring_along.size());
  big_j.reserve(similaring_along.size());
  big_similar.reserve(similaring_along.size());

  // unwrap that nice structure I guess :'(
  while (!similaring_along.empty()) {
    big_similar.push_back(get<0>(similaring_along.top())); // "effing 0-base" is what I should write here :P
    big_i.push_back(get<1>(similaring_along.top()) + 1);
    big_j.push_back(get<2>(similaring_along.top()) + 1);
    similaring_along.pop();
  }


  // make a list object using wrap() for the stdvectors
  return( Rcpp::List::create( 
    Rcpp::Named("big_similar")=wrap( big_similar),
    Rcpp::Named("big_i")=wrap( big_i),
    Rcpp::Named("big_j")=wrap( big_j)
  ));

  END_RCPP
};

// [[Rcpp::export]]
SEXP DUP_paircomps_incomplete_lots(
  RawMatrix geno1, // n_loci, n_samps1,
  RawMatrix geno2, // n_loci, n_samps1,
  bool symmo, // should only be true if geno1==geno2
  double max_diff_ppn, // max allowed ppn discrepant 4way genos for  "true duplicate"
  int limit // exit if ndups hits the limit
) {
  BEGIN_RCPP
    // Avoid scoping woes by doing the include-file right here
#if defined( MVBDEBUG)
#include <r_int_defs_debug.h>
#endif

  int n_samps1;
  int n_samps2;
  int n_loci;
  int n_kept;
  int i;
  int j;
  int g1;
  int g2;
  int j_max;
  int iloc;
  int tot_ndiff;
  int tot_ncomp;
  int is_a_comp;
  int nremloci;

  n_loci = geno1.nrow();
  n_samps1 = geno1.ncol();
  ASSERTO( geno2.nrow() == n_loci);
  n_samps2 = geno2.ncol();

  if( symmo) {
    ASSERTO( n_samps1 == n_samps2);
  };


  // Need to grow the kept-values, so std::vector is apparently better
  std::vector<int> big_ndiff;
  std::vector<int> big_ncomp;
  std::vector<int>  big_i;
  std::vector<int> big_j;

  // check compilation; at one point couldn't find scope
  // now can, but following line doesn't work; prolly doesn't matter

  // All genos should have been recoded s.t. 0 means missing
  n_kept = 0;

  for( i = 0; i < n_samps1; i++){

    j_max = symmo ? i : n_samps2;
    for( j = 0; j < j_max; j++) {
      tot_ncomp = 0;
      tot_ndiff = 0;
      nremloci = n_loci+1;

      for( iloc = 0; iloc < n_loci; iloc++) {
        // Faster version could "unroll" this a bit and only
        // ... do the if-check every 10 loci or so
        // ... and pre-compute a mini-vec of g1's

        g1 = geno1( iloc, i);
        g2 = geno2( iloc, j);
        is_a_comp = (int)( g1 * g2 > 0);
        tot_ncomp += is_a_comp;
        tot_ndiff += is_a_comp * (int)( g1 != g2);

        // Exit if deffo nondup (ie even if all remaining loci match)
        nremloci--;
        if( tot_ndiff > max_diff_ppn * (tot_ncomp + nremloci)) {
      break;
        };
      };

      if( tot_ndiff > max_diff_ppn * tot_ncomp ) {
    continue;
      } else {
        big_ndiff. push_back( tot_ndiff);
        big_ncomp. push_back( tot_ncomp);
        big_i. push_back( i+1); // Effing 0-base...
        big_j. push_back( j+1); // Effing 0-base...
        n_kept += 1;
        if( n_kept >= limit) {
  return( wrap( i)); // how far thru we got. R will detect this isn't a list.
        }; // if limit
      }; // if dup found
    }; // for j
  }; // for i

     // make a list object using wrap() for the stdvectors
  return( Rcpp::List::create(
    Rcpp::Named("big_ndiff")=wrap( big_ndiff),
    Rcpp::Named("big_ncomp")=wrap( big_ncomp),
    Rcpp::Named("big_i")=wrap( big_i),
    Rcpp::Named("big_j")=wrap( big_j)
  ));

  END_RCPP
};
// migrated from MVB 'kinference.cpp', 4/10/18

// [[Rcpp::export]]
SEXP indiv_lglk_geno(
  NumericMatrix lpgeno, // n_loci, n_genotypes
  RawMatrix geno // n_samps, n_loci
) {
  BEGIN_RCPP
    // Avoid scoping woes by doing the include-file right here
#if defined( MVBDEBUG)
#include <../genobasics-vistu/r_int_defs_debug.h>
#endif

  int n_samps;
  int n_loci;
  int i;
  int iloc;
  int this_g;
  double this_lglk;

  n_loci = geno.ncol();
  n_samps = geno.nrow();
  ASSERTO( lpgeno.nrow() == n_loci);

  NumericVector lglk( n_samps);


  // check compilation; at one point couldn't find scope
  // now can, but following line doesn't work; prolly doesn't matter
  // n_kept = XINTEGER(n_PLODs_below)[1];

  for( i = 0; i < n_samps; i++){
    this_lglk = 0;
    for( iloc = 0; iloc < n_loci; iloc++) {
      this_g = geno( i, iloc);
      this_lglk += lpgeno( iloc, geno( i, iloc)-1);
    }; // for iloc
    lglk[ i] = this_lglk;
  }; // for i

     // make a list object using wrap() for the stdvectors
  return lglk;

  END_RCPP
};


// [[Rcpp::export]]
SEXP K_indiv(
  NumericVector tt, // n_t
  RawMatrix geno, // n_loci, n_samps
  NumericVector vec_LOD, // n_loci, n_geno, n_geno
  NumericMatrix Pg // n_loci, n_geno
){
BEGIN_RCPP
    // Avoid scoping woes by doing the include-file right here
#if defined( MVBDEBUG)
#include <../genobasics-vistu/r_int_defs_debug.h>
#endif

  int n_samps;
  int n_geno;
  int n_loci;
  int n_t;
  int i;
  int l;
  int g;
  int i_t;
  int gil;
  double ipdot;
  double temp;
  double plg;
  double this_LOD;

  n_loci = geno.nrow();
  n_samps = geno.ncol();

  IntegerVector LOD_dims= vec_LOD.attr("dim");
  ASSERTO( LOD_dims. size() == 3);
  n_geno = LOD_dims[1];
  ASSERTO( LOD_dims[2] == n_geno);
  ASSERTO( LOD_dims[0] == n_loci);
  arma::cube LOD (vec_LOD.begin(), n_loci, n_geno, n_geno, false);

  ASSERTO( Pg.nrow()==n_loci);
  ASSERTO( Pg.ncol()==n_geno);

  NumericMatrix K( n_t, n_samps);
  NumericMatrix dK( n_t, n_samps);
  NumericMatrix ddK( n_t, n_samps);
  std::fill( K.begin(), K.end(), 0); // shouldn't this be sugarized a la NuMat x(ni,nj)=0 ?
  std::fill( dK.begin(), dK.end(), 0);
  std::fill( ddK.begin(), ddK.end(), 0);
  K = 0; // sugar test

  NumericVector pdot (n_geno);
  NumericVector xpdot (n_geno);
  NumericVector xxpdot (n_geno);

  NumericVector e_CLOD( n_samps);
  NumericVector e2_CLOD( n_samps);
  e_CLOD = 0;
  e2_CLOD = 0;

  for( i= 0; i< n_samps; i++ ) {
    for( l= 0; l< n_loci; l++) {
      gil = geno( i, l);

      for( g= 0; g< n_geno; g++) {
        this_LOD= LOD( l, gil, g);
        plg = Pg( l, g);

        e_CLOD[ i] += plg * this_LOD;
        e2_CLOD[ i] += plg * this_LOD * this_LOD;

        pdot= 0; // vectorized
        xpdot= 0;
        xxpdot= 0;

        for( i_t= 0; i_t< n_t; i_t++) {
          temp = plg * exp( tt[ i_t] * this_LOD);
          pdot[ i_t] += temp;
          temp *= this_LOD;
          xpdot[ i_t] += temp;
          temp *= this_LOD;
          xxpdot[ i_t] += temp;
        }; // i_t
      }; // g

      for( i_t= 0; i_t < n_t; i_t++) {
        K( i_t, i) += log( pdot[ i_t]);
        ipdot = 1/pdot[ i_t];
        temp = xpdot[ i_t] * ipdot;
        dK( i_t, i) += temp;
        ddK( i_t, i) += xxpdot[ i_t] * ipdot - temp*temp;
      }; // i_t
    }; // l
  }; // i

return( Rcpp::List::create( 
  Rcpp::Named("K")=K
));
END_RCPP
};

//    ETT[ it, i, l, g] := exp( tt[ it] * LOD[ l, geno[ i, l], g])
//      log_S[ it, i, l] := log( Pg[ l, g] %[g]% ETT[ it, i, l, g])
//      K[ it, i] := log_S[ it, i, +.]
//SL[ it, l] := PUPLOD[ l, g12] %[g12]% ETT[ it, l, g12]
// dK = rowSums( SL/S);
//SLL[ it, l] := PUPLOD2[ l, g12] %[g12]% ETT[ it, l, g12]
// rowSums( (SLL/S-sqr( SL/S)))

