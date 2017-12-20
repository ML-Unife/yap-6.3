#include <cassert>

#include "GenericFactor.h"
#include "ProbFormula.h"
#include "Indexer.h"


namespace Horus {

template <typename T> const T&
GenericFactor<T>::argument (size_t idx) const
{
  assert (idx < args_.size());
  return args_[idx];
}



template <typename T> T&
GenericFactor<T>::argument (size_t idx)
{
  assert (idx < args_.size());
  return args_[idx];
}



template <typename T> unsigned
GenericFactor<T>::range (size_t idx) const
{
  assert (idx < ranges_.size());
  return ranges_[idx];
}



template <typename T> bool
GenericFactor<T>::contains (const T& arg) const
{
  return Util::contains (args_, arg);
}



template <typename T> bool
GenericFactor<T>::contains (const std::vector<T>& args) const
{
  for (size_t i = 0; i < args.size(); i++) {
    if (contains (args[i]) == false) {
      return false;
    }
  }
  return true;
}



template <typename T> void
GenericFactor<T>::setParams (const Params& newParams)
{
  params_ = newParams;
  assert (params_.size() == Util::sizeExpected (ranges_));
}



template <typename T> double
GenericFactor<T>::operator[] (size_t idx) const
{
  assert (idx < params_.size());
  return params_[idx];
}



template <typename T> double&
GenericFactor<T>::operator[] (size_t idx)
{
  assert (idx < params_.size());
  return params_[idx];
}



template <typename T> GenericFactor<T>&
GenericFactor<T>::multiplyHet (const GenericFactor<T>& g, size_t idxConv1,
size_t idxConv2)
{
  //std::cout <<"het mult "<<idxConv1<< " " <<idxConv2<<"\n";
  unsigned range_prod = 1;
//  std::vector<T> g_args_or = g.arguments();
//  Ranges g_ranges_or  = g.ranges();
//  std::vector<T> args1or=args_;
//  Ranges ranges1or=ranges_;
//  Params params1or=params_;
   std::vector<T> g_args = g.arguments();
  Ranges g_ranges  = g.ranges();
  Params g_params  = g.params();
  const std::vector<bool> g_conv = g.convergent();
  std::vector<T> args1=args_;
  Ranges ranges1=ranges_;
  Params params1=params_;
  
  for (size_t i = 0; i < g_args.size(); i++) {
    size_t idx = indexOf (g_args[i]);
  //  std::cout<< idx<<" idx\n";
    if (idx == args_.size()) {
      range_prod *= g_ranges[i];
      args_.push_back (g_args[i]);
      ranges_.push_back (g_ranges[i]);
      convergent_.push_back(g_conv[i]);
    } else {
      convergent_[idx]=convergent_[idx]||g_conv[i];
    //std::cout<< "share\n";
    }
  }
  std::vector<T> args_mult=args_;
  Ranges ranges_mult=ranges_;
  T argconv=args_mult[idxConv1];
  args_mult.erase(args_mult.begin()+idxConv1);
  args_mult.push_back(argconv);
  size_t rangeconv = ranges_mult[idxConv1];
  ranges_mult.erase(ranges_mult.begin()+idxConv1);
  ranges_mult.push_back(rangeconv);
/*  for (size_t i = 0; i < args_.size(); i++) {
    if (i != 0) std::cout << ", " ;
    std::cout << args_[i].print<<"ciao";
  }
  */

    //T arg1conv=args1[idxConv1];
    //args1.erase(args1.begin()+idxConv1);
    //args1.push_back(arg1conv);
    //size_t range1conv=ranges1[idxConv1];
    //ranges1.erase(ranges1.begin()+idxConv1);
    //ranges1.push_back(range1conv);
    


//    T arg2conv=g_args[idxConv2];
 //   g_args.erase(g_args.begin()+idxConv2);
 //   g_args.push_back(arg2conv);
//    size_t range2conv=g_ranges[idxConv2];
//    g_ranges.erase(g_ranges.begin()+idxConv2);
//    g_ranges.push_back(range2conv);
    
    //size_t size=params_.size()*range_prod;
    //Params p1 (size,0);
    //Params p2 (size,0);
    //Params p (size,0);
// std::cout<< Horus::Parfactor::Parfactor (args1, p1, tuples, distId_, convergent_);

    extend (range_prod);
    MapIndexer indexer ( args_mult, ranges_mult,args_, ranges_);
    MapIndexer indexer1 (args_mult, ranges_mult,args1, ranges1);
    MapIndexer indexer2 (args_mult, ranges_mult,g_args, g_ranges);
    if (Globals::logDomain) {
      for (; indexer.valid();  ) {
        
      }
    } else {


      for (; indexer.valid(); ) {
         double c1f=params1[indexer1];
         double c2f=g_params[indexer2];
         params_[indexer]= c1f*c2f;
         ++indexer;
	 ++indexer1;
	 ++indexer2;
         double c1t=params1[indexer1];
         double c2t=g_params[indexer2];
//std::cout<<indexer <<"\n";
 	params_[indexer] =c1f*c2t+c1t*c2f+c1t*c2t;
	++indexer;
        ++indexer1;
         ++indexer2;
 //	std::cout<< c1f<<" "<<c1t<<" "<<c2f<<" "<<c2t<<"\n"<<p<<"\n";
	 }
}
 /*
    MapIndexer indexerp1 ( args_mult, ranges_mult,args_, ranges_);
    MapIndexer indexer1 (args_mult, ranges_mult,args1, ranges1);
    if (Globals::logDomain) {
      for (; indexerp1.valid();  ) {
        
      }
    } else {
      for (; indexerp1.valid(); ) {
      std::cout<<indexerp1<<" indexer\n"<<indexer1<<" indexer1\n";
         p1[indexerp1]= params1[indexer1];
	 ++indexer1;
	 ++indexerp1;
      std::cout<<indexerp1<<" indexer\n"<<indexer1<<" indexer1\n";
 	 p1[indexerp1]= params1[indexer1];
	 ++indexer1;
	 ++indexerp1;
	 }
    std::cout << p1<<"p1\n";
    MapIndexer indexerp2(args_mult, ranges_mult, args_, ranges_);
    MapIndexer indexer2 (args_mult, ranges_mult, g_args, g_ranges);

	       for (; indexerp2.valid(); ) {
      std::cout<<indexerp2<<" indexer\n"<<indexer2<<" indexer2\n";
//      std::cout<<indexer<<" indexer\n"<<indexer1<<" indexer1\n"<<indexer2<<"ind2\n";
         p2[indexerp2]= g_params[indexer2];
         ++indexer2;
         ++indexerp2;
      std::cout<<indexerp2<<" indexer\n"<<indexer2<<" indexer2\n";
//      std::cout<<indexer<<" indexer\n"<<indexer1<<" indexer1\n"<<indexer2<<"ind2\n";
         p2[indexerp2]= g_params[indexer2];
         ++indexer2;
         ++indexerp2;
         }
    std::cout << p2<<"p1\n";
    MapIndexer indexer(args_mult, ranges_mult, args_, ranges_);
    ;

	       for (; indexer.valid(); ) {
      std::cout<<indexer<<" indexer\n";
         double c1f=p1[indexer];
         double c2f=p2[indexer];
         p[indexer]= c1f*c2f;
         ++indexer;
         double c1t=p1[indexer];
         double c2t=p2[indexer];
std::cout<<indexer <<"\n";
 	p[indexer] =c1f*c2t+c1t*c2f+c1t*c2t;
	++indexer;
 	std::cout<< c1f<<" "<<c1t<<" "<<c2f<<" "<<c2t<<"\n"<<p<<"\n";
      }
    }*/
//    params_=p;
  return *this;
}


template <typename T> GenericFactor<T>&
GenericFactor<T>::multiplyHet2 (const GenericFactor<T>& g, std::vector<size_t> idxConv1,
std::vector<size_t> idxConv2)
{
  //std::cout <<"het mult 2 "<<idxConv1<< " " <<idxConv2<<"\n";
  unsigned range_prod = 1;
//  std::vector<T> g_args_or = g.arguments();
//  Ranges g_ranges_or  = g.ranges();
//  std::vector<T> args1or=args_;
//  Ranges ranges1or=ranges_;
//  Params params1or=params_;
  
   std::vector<T> g_args = g.arguments();
  Ranges g_ranges  = g.ranges();
  Params g_params  = g.params();
  const std::vector<bool> g_conv = g.convergent();
  
  std::vector<T> args1=args_;
  Ranges ranges1=ranges_;
  Params params1=params_;
  
  size_t sizeConv1 = 0;
  size_t sizeConv2 = 0;
 
  //std::cout<<convergent_<<"\n";
  for (size_t i = 0; i < convergent_.size(); i++) {
     if (convergent_[i]) sizeConv1++;
  }
  //std::cout<<g_conv<<"\n";
  for (size_t i = 0; i < g_conv.size(); i++) {
     if (g_conv[i]) sizeConv2++;
  }
 
  for (size_t i = 0; i < g_args.size(); i++) {
    size_t idx = indexOf (g_args[i]);
  //  std::cout<< idx<<" idx\n";
    if (idx == args_.size()) {
      range_prod *= g_ranges[i];
      args_.push_back (g_args[i]);
      ranges_.push_back (g_ranges[i]);
      convergent_.push_back(g_conv[i]);
    } else {
      convergent_[idx]=convergent_[idx]||g_conv[i];
      //share = true;
      //std::cout<< "share\n";
    }
  }
  
  std::vector<T> args_mult=args_;
  Ranges ranges_mult=ranges_;
  for (size_t idx = 0; idx < idxConv1.size(); idx++) {
        T argconv=args_mult[idxConv1[idx]];
        args_mult.erase(args_mult.begin()+idxConv1[idx]-idx);
        args_mult.push_back(argconv);
        size_t rangeconv = ranges_mult[idxConv1[idx]];
        ranges_mult.erase(ranges_mult.begin()+idxConv1[idx]-idx);
        ranges_mult.push_back(rangeconv);
  }
  
  
  
    extend (range_prod);
  
    MapIndexer indexer ( args_mult, ranges_mult,args_, ranges_);
    MapIndexer indexer1 (args_mult, ranges_mult,args1, ranges1);
    MapIndexer indexer2 (args_mult, ranges_mult,g_args, g_ranges);
    
    std::vector<bool> curr_label ((sizeConv1 > sizeConv2) ? 
                                                sizeConv1 : sizeConv2, false);
    
    
    size_t reset_position_1 = 0;
    size_t reset_position_2 = 0;
    size_t n_res = 0;
    size_t lim1 = std::pow(2,sizeConv1), lim2 = std::pow(2,sizeConv2), lim = std::pow(2,curr_label.size());
    
    //std::cout<<"lim1: "<<lim1<<" lim2: "<<lim2<<" lim: "<<lim<<"\n";
    for (;indexer.valid(); n_res++) {
        
        std::vector<bool> label1 (sizeConv1, false);
        std::vector<bool> label2 (sizeConv2, false);
        
        std::vector<std::vector<bool>> labelsEl1, labelsEl2;
        std::vector<double> el1, el2;
        
        //std::cout<<"n_res "<<n_res<<"\n";
        if (n_res == lim) {
            //std::cout<<"reset\n";
            n_res = 0;
            reset_position_1 += lim1;
            reset_position_2 += lim2;
            curr_label = resetLabel(curr_label);
        }
        
        indexer1.reset();
        for (size_t i1 = 0; i1 < reset_position_1; ++indexer1,i1++);
        indexer2.reset();
        for (size_t i2 = 0; i2 < reset_position_2; ++indexer2,i2++);
        
        
        //prendo i valori di this
        for (size_t i = 0; i < lim1 && curr_label[0] >= label1[0]; i++) {
            if (lessLabels(curr_label,label1)) {
                labelsEl1.push_back(label1);
                el1.push_back(params1[indexer1]);
                ++indexer1;
                label1=incrementLabel(label1);
            } else {
                ++indexer1;
            }
        }
        //std::cout<<labelsEl1<<"\n"<<el1<<"\n";
        
        //prendo i valori di g
        for (size_t i = 0; i < lim2 && curr_label[0] >= label2[0]; i++) {
            if (lessLabels(curr_label,label2)) {
                labelsEl2.push_back(label2);
                el2.push_back(g_params[indexer2]);
                ++indexer2;
                label2=incrementLabel(label2);
            } else {
                ++indexer2;
            }
        }
        //std::cout<<labelsEl2<<"\n"<<el2<<"\n";
        
        //std::cout<<"mult for "<<curr_label<<"\n";
        params_[indexer] = 0;
        for (size_t i = 0; i < el1.size(); i++) {
            for (size_t j = 0; j < el2.size(); j++) {
                if (orLabels(labelsEl1[i], labelsEl2[j], curr_label)) {
                    //std::cout<<labelsEl1[i]<<" + "<<labelsEl2[j]<<"\n";
                    params_[indexer] += el1[i]*el2[j];
                }
            }
        }
        
        ++indexer;
        curr_label = incrementLabel(curr_label);
        
    }
    
    
  return *this;
}

template <typename T> GenericFactor<T>&
GenericFactor<T>::multiply (const GenericFactor<T>& g)
{
  if (args_ == g.arguments()) {
    // optimization
    Globals::logDomain
      ? params_ += g.params()
      : params_ *= g.params();
    return *this;
  }
  unsigned range_prod = 1;
  bool share_arguments = false;
  const std::vector<T>& g_args = g.arguments();
  const Ranges& g_ranges  = g.ranges();
  const Params& g_params  = g.params();
  const std::vector<bool> g_conv = g.convergent();
  for (size_t i = 0; i < g_args.size(); i++) {
    size_t idx = indexOf (g_args[i]);
  //  std::cout<< idx<<" idx\n";
    if (idx == args_.size()) {
      range_prod *= g_ranges[i];
      args_.push_back (g_args[i]);
      ranges_.push_back (g_ranges[i]);
      convergent_.push_back(g_conv[i]);
    } else {
      convergent_[idx]=convergent_[idx]||g_conv[i];
    //std::cout<< "share\n";
      share_arguments = true;
    }
  }
  if (share_arguments == false) {
    // optimization
    cartesianProduct (g_params.begin(), g_params.end());
  } else {
    extend (range_prod);
    Params::iterator it = params_.begin();
    MapIndexer indexer (args_, ranges_, g_args, g_ranges);
    if (Globals::logDomain) {
      for (; indexer.valid(); ++it, ++indexer) {
        *it += g_params[indexer];
      }
    } else {
      for (; indexer.valid(); ++it, ++indexer) {
        *it *= g_params[indexer];
      }
    }
  }
  return *this;
}



template <typename T> void
GenericFactor<T>::sumOutIndexNOV (size_t idx, size_t idxConv, size_t count)
{
  assert (idx < args_.size());
  assert (args_.size() > 1);
  size_t idxConvRes=idxConv;
  size_t size_res = params_.size() / ranges_[idx];
  std::vector<T>  args_no_conv=args_;
  std::vector<T>  args_res=args_;
  //std::cout<< "het sum idx "<< idx<<" idxcon "<< idxConv<<" count "<< count<<" args " << args_.size() <<"\n";
  Ranges ranges_no_conv=ranges_;
  Ranges ranges_res=ranges_;
  std::vector<bool> conv_res=convergent_;
  T arg=args_[idx];
  T argConv=args_[idxConv];
  //std::cout<<" arg " << arg <<" argConv " << argConv<<"\n";
  std::vector<T>  args_or=args_;
  args_or.erase(args_or.begin()+idx);
  //std::cout<<"or\n";
  //std::cout<<args_or<<"\n";
  Ranges ranges_or=ranges_;
  ranges_or.erase(ranges_or.begin()+idx);
  //std::cout<<ranges_or<<"\n";
  if (idxConv>idx)
  {
    idxConvRes=idxConv-1;
    args_no_conv.erase (args_no_conv.begin() + idx);
    args_no_conv.erase (args_no_conv.begin() + idxConvRes);
    ranges_no_conv.erase (ranges_no_conv.begin() + idx);
    ranges_no_conv.erase (ranges_no_conv.begin() + idxConvRes);
   }
  else {
//  if (idxConv < idx)
//  { 
    args_no_conv.erase (args_no_conv.begin() + idxConv);
    args_no_conv.erase (args_no_conv.begin() + idx-1);
    ranges_no_conv.erase (ranges_no_conv.begin() + idxConv);
    ranges_no_conv.erase (ranges_no_conv.begin() + idx-1);
//  } else {
//    args_no_conv.erase (args_no_conv.begin() + idx);
//    ranges_no_conv.erase (ranges_no_conv.begin() + idx);
//    ranges_no_conv.erase (ranges_no_conv.begin() + idxConv);
//    args_no_conv.erase (args_no_conv.begin() + idxConv);
//  }
  }
  //std::cout<<"no_conv\n";
  //std::cout<<args_no_conv<<"\n";
  //std::cout<<ranges_no_conv<<"\n";
  args_res.erase(args_res.begin()+idx);
  ranges_res.erase(ranges_res.begin()+idx);
  conv_res.erase(conv_res.begin()+idx);
  //std::cout<<"res\n";
  //std::cout<<args_res<<"\n";
  //std::cout<<ranges_res<<"\n";
  //std::cout<<conv_res<<"\n";
  //std::cout <<"par "<<params_;
/*  args.push_back(arg1);
  args.push_back(arg);*/
 
//  reorderArguments(args);
  //std::cout<<args<<"args\n";
// std::cout <<"par "<<params_;
  
  Params ps_res (size_res, LogAware::addIdenty());
 // Params::iterator np  = newps.begin();


  std::vector<T>  args_reordered=args_no_conv;
  Ranges ranges_reordered=ranges_no_conv;
  args_reordered.push_back(argConv);
  args_reordered.push_back(arg);
  ranges_reordered.push_back(ranges_[idxConv]);
  ranges_reordered.push_back(ranges_[idx]);
 
// std::cout <<"\n\n"<< args1.size() << "\n\n";
  MapIndexer indexer_reordering (args_reordered,ranges_reordered, args_,ranges_);
MapIndexer indexer_res (args_res,ranges_res, args_or,ranges_or);

//std::cout<<"res\n";
//std::cout<<ps_res<<"\n";


    for (; indexer_reordering.valid();) {

//    std::cout<<last-first<<" l-f\n"<< indexer<<" ind\n";
//      std::cout<< newps[indexer] <<" newps "<<*first <<" first\n"; 
      double a= params_[indexer_reordering];
      ++indexer_reordering;
      double b= params_[indexer_reordering];
      ++indexer_reordering;
      double c= params_[indexer_reordering];
      ++indexer_reordering;
      double d= params_[indexer_reordering];
      ++indexer_reordering;
      //std::cout<<"a "<<a<<" b "<<b<<" c "<<c<<" d "<<d<<"\n";
      
      double ab= LogAware::pow(a+b,(unsigned)count);
      ps_res[indexer_res]=ab;
      ++indexer_res;
      ps_res[indexer_res]=LogAware::pow(a+b+c+d,(unsigned)count)-ab;
      ++indexer_res;
      //std::cout<<"res\n";
        //std::cout<<ps_res<<"\n";
     // std::cout<<params_<<"\n";
      //std::cout<<newps<<"\n\n";
    }
  
   ranges_=ranges_res;
//  std::cout<<"fine-1\n";
  params_ = ps_res;
  args_=args_res;
  convergent_= conv_res;
  //std::cout<<"fine0 size arg " << args_.size() <<"\n";
  if (args_.size()==1 &&  heterogeneous_)
  {
    //std::cout << "only one conv var\n";
    heterogeneous_ = false;
    //std::cout << "beg "<< *convergent_.begin();
    //std::cout<< "convergent "<< convergent_<<"here\n";
    //std::cout.flush();
    convergent_[0]=false;
  } 

  
  //std::cout<<"fine\n";
}

template <typename T> void
GenericFactor<T>::sumOutIndexNOV2 (size_t idx, std::vector<size_t> idxConv, size_t count)
{
  assert (idx < args_.size());
  assert (args_.size() > 1);
  
  size_t size_res = params_.size() / ranges_[idx];
  
  // args vectors
  std::vector<T>  args_no_conv;
  std::vector<T>  args_conv;
  std::vector<T>  args_res=args_;
  
  //std::cout<< "het sum 2 idx "<< idx<<" idxcon "<< idxConv<<" count "<< count<<" args " << args_.size() <<"\n";
  
  //ranges vector
  Ranges ranges_no_conv;
  Ranges ranges_conv;
  Ranges ranges_res=ranges_;
  
  std::vector<bool> conv_res=convergent_;
  
  T arg=args_[idx];

  
  //std::cout<<" arg " << arg <<" argConv " << argConv<<"\n";
  // or vectors
  std::vector<T>  args_or=args_;
  args_or.erase(args_or.begin()+idx);
  
  Ranges ranges_or=ranges_;
  ranges_or.erase(ranges_or.begin()+idx);
  
  //print==
  //std::cout<<"or\n";
  //std::cout<<ranges_or<<"\n";
//  for (size_t i = 0; i < args_or.size(); i++) {
//    if (i != 0) std::cout << ", " ;
//    std::cout << args_or[i];
//  }
  //==
  
  for (size_t i = 0; i < args_.size(); i++) {
      if (convergent_[i] && i != idx) {
          //std::cout<<"conv "<<args_conv.size()<<"\n";
        args_conv.push_back(args_[i]);
        ranges_conv.push_back(ranges_[i]);
      } else if (!convergent_[i] && i != idx) {
          //std::cout<<"no-conv "<<args_no_conv.size()<<"\n";
        args_no_conv.push_back(args_[i]);
        ranges_no_conv.push_back(ranges_[i]);
      }
  }
  
  
  
  //res vectors
  
  //print===
  //std::cout<<"no_conv"<<"\n";
  //std::cout<<ranges_no_conv<<"\n";
//  for (size_t i = 0; i < args_no_conv.size(); i++) {
//    if (i != 0) std::cout << ", " ;
//    std::cout << args_no_conv[i];
//  }
  //==
  
  // res vector
  args_res.erase(args_res.begin()+idx);
  ranges_res.erase(ranges_res.begin()+idx);
  conv_res.erase(conv_res.begin()+idx);
  
  //print==
  //std::cout<<"res\n";
  //std::cout<<ranges_res<<"\n";
  //std::cout<<conv_res<<"\n";
//  for (size_t i = 0; i < args_res.size(); i++) {
//    if (i != 0) std::cout << ", " ;
//    std::cout<< args_res[i];
//  }
  //==
  
  Params ps_res (size_res, LogAware::addIdenty());
 // Params::iterator np  = newps.begin();


  std::vector<T>  args_reordered;
  args_reordered.reserve(args_.size());//=args_no_conv;
  Ranges ranges_reordered;
  ranges_reordered.reserve(ranges_.size()); //=ranges_conv;
  args_reordered.insert(args_reordered.end(), args_no_conv.begin(), args_no_conv.end());
  //std::cout<<"reord "<<args_reordered.size()<<"\n";
  args_reordered.insert(args_reordered.end(), args_conv.begin(), args_conv.end());
  //std::cout<<"reord "<<args_reordered.size()<<"\n";
  args_reordered.push_back(arg);
  //std::cout<<"reord "<<args_reordered.size()<<"\n";
  ranges_reordered.insert(ranges_reordered.end(), ranges_no_conv.begin(), ranges_no_conv.end());
  ranges_reordered.insert(ranges_reordered.end(), ranges_conv.begin(), ranges_conv.end());
  ranges_reordered.push_back(ranges_[idx]);
 
  //print==
  //std::cout<<"reordered\n";
  //std::cout<<args_res<<"\n";
  //std::cout<<ranges_reordered<<"\n";
//    for (size_t i = 0; i < args_reordered.size(); i++) {
//    if (i != 0) std::cout << ", " ;
//    std::cout << args_reordered[i];
//  }
  //==
  
  // std::cout <<"\n\n"<< args1.size() << "\n\n";
  MapIndexer indexer_reordering (args_reordered,ranges_reordered, args_,ranges_);
  MapIndexer indexer_res (args_res,ranges_res, args_or,ranges_or);
  MapIndexer indexer_sottr_res (args_res,ranges_res, args_or,ranges_or);

  //std::cout<<"res\n";
  //std::cout<<ps_res<<"\n";

  size_t n_elements= 1;
  std::vector<bool> curr_label (args_conv.size(), false); //args_conv
  std::vector<bool> label (args_conv.size(), false);
  //std::cout<<"cl "<<curr_label<<" l "<<label<<"\n";
  
  // for each position in res
  size_t reset_position_reordering = 0;
  size_t reset_position_sottr = 0;
  for (size_t res_idx = 0; res_idx < ps_res.size(); res_idx++) {
     
            double res = 0;
            indexer_reordering.reset();
            //std::cout<<"aumento index_reordering\n";
            for (size_t irpr = 0; irpr < reset_position_reordering; ++indexer_reordering,irpr++);
            //std::cout<<"fine aumento index_reordering\n";
            //check n_elements*2
            //std::cout<<"n_el "<<n_elements<<"\n";
            for (size_t j = 0; j < n_elements; j++) {
                if (indexer_reordering.valid()) {

                    //std::cout<<"cl "<<curr_label<<" l "<<label<<"\n";
                    if (lessLabels(curr_label, label)) {
                        //std::cout<<"less\n";

                      res += params_[indexer_reordering];
                      ++indexer_reordering;

                      res += params_[indexer_reordering];
                      ++indexer_reordering;

                    } else {
                        ++indexer_reordering;
                        ++indexer_reordering;
                    }
                    label=incrementLabel(label);
                }
            }
            n_elements++;

            ps_res[indexer_res] = LogAware::pow(res,(unsigned) count);
            //std::cout<<"\n"<<ps_res<<"\n";
            // sottrai minori
            label=resetLabel(label);
            indexer_sottr_res.reset();//[];
            //std::cout<<"aumento index_sottr\n";
            for (size_t irps = 0; irps < reset_position_sottr; ++indexer_sottr_res,irps++);
            //std::cout<<"aumento index_sottr\n";

            //inserimento nel vettore risultati
            for (size_t j = reset_position_sottr; j < res_idx; j ++) {
                //std::cout<<"i "<<res_idx<<" j "<<j<<"\n";
                if (lessLabels(curr_label, label)) {
                    //std::cout<<"cl "<<curr_label<<" l "<<label<<"\n";
                    ps_res[indexer_res] -= ps_res[indexer_sottr_res];
                    //std::cout<<" l "<<label<<"\n";
                    label=incrementLabel(label);
                    //std::cout<<" l "<<label<<"\n";
                    ++indexer_sottr_res;
                }
            }
            //std::cout<<"\n"<<ps_res<<"\n";
            

            label=resetLabel(label);
            curr_label=incrementLabel(curr_label);
            ++indexer_res;
            
            //std::cout<<"n_elements " <<n_elements<<"\n";
            
            float lim = std::pow(2,args_conv.size());
            if (n_elements > lim) {
                //std::cout<<"reset\n";
                //std::cout<<"nci "<<res_idx<<" n_el "<<n_elements;
                
                reset_position_reordering += lim * 2;
                reset_position_sottr = res_idx + 1;
                //resettare label
                n_elements = 1;
                curr_label = resetLabel(curr_label);
                //controllare
            }
      
      
  }
  
  //std::cout<<"\n"<<ps_res<<"\n";
//  indexer_reordering.reset();
//  indexer_res.reset();
  
//    for (; indexer_reordering.valid();) {
//
////    std::cout<<last-first<<" l-f\n"<< indexer<<" ind\n";
////      std::cout<< newps[indexer] <<" newps "<<*first <<" first\n"; 
//      double a= params_[indexer_reordering];
//      std::cout<<"a";
//      ++indexer_reordering;
//      double b= params_[indexer_reordering];
//      std::cout<<"b";
//      ++indexer_reordering;
//      double c= params_[indexer_reordering];
//      std::cout<<"c";
//      ++indexer_reordering;
//      double d= params_[indexer_reordering];
//      std::cout<<"d";
//      ++indexer_reordering;
//      std::cout<<"a "<<a<<" b "<<b<<" c "<<c<<" d "<<d<<"\n";
//      
//      double ab= LogAware::pow(a+b,count);
//      ps_res[indexer_res]=ab;
//      ++indexer_res;
//      ps_res[indexer_res]=LogAware::pow(a+b+c+d,count)-ab;
//      ++indexer_res;
//      std::cout<<"res\n";
//        std::cout<<ps_res<<"\n";
//     // std::cout<<params_<<"\n";
//      //std::cout<<newps<<"\n\n";
//    }
  
  
   ranges_=ranges_res;
//  std::cout<<"fine-1\n";
  params_ = ps_res;
  args_=args_res;
  convergent_= conv_res;
  //std::cout<<"fine0 size arg " << args_.size() <<"\n";
  if (args_.size()==1 &&  heterogeneous_)
  {
    //std::cout << "only one conv var\n";
    heterogeneous_ = false;
    //std::cout << "beg "<< *convergent_.begin();
    //std::cout<< "convergent "<< convergent_<<"here\n";
    //std::cout.flush();
    convergent_[0]=false;
  } 

  
  //std::cout<<"fine\n";
}

template <typename T> std::vector<bool>
GenericFactor<T>::incrementLabel(std::vector<bool> label) {
    for (int i = label.size() - 1; i >= 0 ; i--) {
        if (label[i])
            label[i] = false;
        else {
            label[i] = true;
            break;
        }
    }
    return label;
}

template <typename T> std::vector<bool>
GenericFactor<T>::resetLabel(std::vector<bool> label) {
    for (int i = label.size() - 1; i >= 0 ; i--) {
        label[i] = false;
    }
    return label;
}

// B < A
template <typename T> bool
GenericFactor<T>::lessLabels(std::vector<bool> labelA, std::vector<bool> labelB) {
    size_t idx_comp;
    if (labelA.size() >= labelB.size()) {
        for (idx_comp = 0; ((idx_comp < labelB.size()) && 
                (labelA[idx_comp] >= labelB[idx_comp])); idx_comp++);
        return (idx_comp >= labelB.size());
    } else {
        return false;
    }
}

template <typename T> bool
GenericFactor<T>::orLabels(std::vector<bool> labelA, std::vector<bool> labelB, std::vector<bool> labelRes) {
    size_t idx_comp;
    size_t sizeA = labelA.size(), sizeB = labelB.size();
    for (idx_comp = 0; idx_comp < labelRes.size(); idx_comp++) {
        if (idx_comp < sizeA && idx_comp < sizeB) {
                if (!((labelA[idx_comp] | labelB[idx_comp]) == labelRes[idx_comp]))
                    return false;
        } else if (idx_comp < sizeA) {
            if (!(labelA[idx_comp] < labelRes[idx_comp])) {
                return false;
            }
        } else {
            if (!(labelB[idx_comp] < labelRes[idx_comp])){
                return false;
            }
        }
    }
    return true;
}

template <typename T> void
GenericFactor<T>::sumOutIndex (size_t idx)
{
  assert (idx < args_.size());
  assert (args_.size() > 1);
  

  size_t new_size = params_.size() / ranges_[idx];
  Params newps (new_size, LogAware::addIdenty());
  Params::const_iterator first = params_.begin();
  Params::const_iterator last  = params_.end();
  MapIndexer indexer (ranges_, idx);
  if (Globals::logDomain) {
    for (; first != last; ++indexer) {
      newps[indexer] = Util::logSum (newps[indexer], *first++);
    }
  } else {
    for (; first != last; ++indexer) {
      newps[indexer] += *first++;
    }
  }
  params_ = newps;
  args_.erase (args_.begin() + idx);
  ranges_.erase (ranges_.begin() + idx);
  convergent_.erase(convergent_.begin()+idx);
//  heterogeneous_ = false;
//  for (size_t i = 0; i < convergent_.size(); i++)
//      if (convergent_[i] == true) {
//          heterogeneous_ = true;
//          break;
//      }
}

template <typename T> void
GenericFactor<T>::setHeterogeneous()
{
  heterogeneous_ = false;
  for (size_t i = 0; i < convergent_.size(); i++)
      if (convergent_[i] == true) {
          heterogeneous_ = true;
          break;
      }
}

template <typename T> void
GenericFactor<T>::absorveEvidence (const T& arg, unsigned obsIdx)
{
  size_t idx = indexOf (arg);
  assert (idx != args_.size());
  assert (obsIdx < ranges_[idx]);
  Params newps;
  newps.reserve (params_.size() / ranges_[idx]);
  Indexer indexer (ranges_);
  for (unsigned i = 0; i < obsIdx; ++i) {
    indexer.incrementDimension (idx);
  }
  while (indexer.valid()) {
    newps.push_back (params_[indexer]);
    indexer.incrementExceptDimension (idx);
  }
  params_ = newps;
  args_.erase (args_.begin() + idx);
  ranges_.erase (ranges_.begin() + idx);
}



template <typename T> void
GenericFactor<T>::reorderArguments (const std::vector<T>& new_args)
{
//  std::cout <<"reorder\n";
  assert (new_args.size() == args_.size());
  if (new_args == args_) {
    return; // already on the desired order
  }
  Ranges new_ranges;
  std::vector<bool> new_conv;
  for (size_t i = 0; i < new_args.size(); i++) {
    size_t idx = indexOf (new_args[i]);
    assert (idx != args_.size());
    new_ranges.push_back (ranges_[idx]);
    new_conv.push_back(convergent_[idx]);
   std::cout << i<< "i "<<idx<<" idx\n";

  }
  /*for (size_t i=0;i< convergent_.size();i++)
  {
    size_t idx = Util::indexOf (new_args,args_[convergent_[i]]);
  //  std::cout << idx << " idx "<<convergent_[i]<<"\n";
    new_conv.push_back(idx);
  }*/

  Params newps;
  newps.reserve (params_.size());
  MapIndexer indexer (new_args, new_ranges, args_, ranges_);
  for (; indexer.valid(); ++indexer) {
    newps.push_back (params_[indexer]);
  }
  params_ = newps;
  args_   = new_args;
  ranges_ = new_ranges;
  convergent_= new_conv;
}



template <typename T> void
GenericFactor<T>::extend (unsigned range_prod)
{
  Params backup = params_;
  params_.clear();
  params_.reserve (backup.size() * range_prod);
  Params::const_iterator first = backup.begin();
  Params::const_iterator last  = backup.end();
  for (; first != last; ++first) {
    for (unsigned reps = 0; reps < range_prod; ++reps) {
      params_.push_back (*first);
    }
  }
}



template <typename T> void
GenericFactor<T>::cartesianProduct (
  Params::const_iterator first2,
  Params::const_iterator last2)
{
  Params backup = params_;
  params_.clear();
  params_.reserve (params_.size() * (last2 - first2));
  Params::const_iterator first1 = backup.begin();
  Params::const_iterator last1  = backup.end();
  Params::const_iterator tmp;
  if (Globals::logDomain) {
    for (; first1 != last1; ++first1) {
      for (tmp = first2; tmp != last2; ++tmp) {
        params_.push_back ((*first1) + (*tmp));
      }
    }
  } else {
    for (; first1 != last1; ++first1) {
      for (tmp = first2; tmp != last2; ++tmp) {
        params_.push_back ((*first1) * (*tmp));
      }
    }
  }
}

template <typename T>
const bool GenericFactor<T>::isDeputy()
{
  bool dep=false;
  for (size_t i=0;i<convergent_.size();i++)
     dep=dep || convergent_[i];
  return dep && !heterogeneous_;
}


template class GenericFactor<VarId>;
template class GenericFactor<ProbFormula>;

}  // namespace Horus

