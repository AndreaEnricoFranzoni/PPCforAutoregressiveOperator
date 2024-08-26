#ifndef KO_DEFAULT_PARAMS_HPP
#define KO_DEFAULT_PARAMS_HPP


namespace DEF_PARAMS_PPC
{

enum CV_type
{
  NoCV = 0,
  CV   = 1,
};

//type for which strategy is used to remove nans
enum MA_type
{
  EMA = 0,      //replacing nans with exponential moving average
  WMA = 1,      //replacing nans with weighted moving average
  SMA = 2,      //replacing nans with simple moving average
  MR  = 3,      //replacing nans with mean (easily changes the mean of the distribution)
  ZR  = 4,      //replacing nans with 0s (easily changes the sd of the distribution)
};


} //end namespace DEF_PARAMS_PPC

#endif /*KO_DEFAULT_PARAMS_HPP*/