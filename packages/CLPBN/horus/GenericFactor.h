#ifndef YAP_PACKAGES_CLPBN_HORUS_GENERICFACTOR_H_
#define YAP_PACKAGES_CLPBN_HORUS_GENERICFACTOR_H_

#include <vector>

#include "Util.h"


namespace Horus {

template <typename T>
class GenericFactor {
  public:
    const std::vector<T>& arguments() const { return args_; }

    std::vector<T>& arguments() { return args_; }

    const Ranges& ranges() const { return ranges_; }

    const Params& params() const { return params_; }

    Params& params() { return params_; }

    size_t nrArguments() const { return args_.size(); }

    size_t size() const { return params_.size(); }

    unsigned distId() const { return distId_; }

    void setDistId (unsigned id) { distId_ = id; }

    void normalize() { LogAware::normalize (params_); }

    size_t indexOf (const T& t) const { return Util::indexOf (args_, t); }

    const T& argument (size_t idx) const;

    T& argument (size_t idx);

    unsigned range (size_t idx) const;

    bool contains (const T& arg) const;

    bool contains (const std::vector<T>& args) const;

    void setParams (const Params& newParams);

    double operator[] (size_t idx) const;

    double& operator[] (size_t idx);

    GenericFactor<T>& multiply (const GenericFactor<T>& g);

     GenericFactor<T>& multiplyHet (const GenericFactor<T>& g, size_t idxConv1,
size_t idxConv2);
     
     GenericFactor<T>& multiplyHet2 (const GenericFactor<T>& g, std::vector<size_t> idxConv1,
std::vector<size_t> idxConv2);

    void sumOutIndex (size_t idx);

    void sumOutIndexNOV (size_t idx, size_t idxConv, size_t count);
    
    void sumOutIndexNOV2 (size_t idx, std::vector<size_t> idxConv, size_t count);

    void absorveEvidence (const T& arg, unsigned obsIdx);

    void reorderArguments (const std::vector<T>& new_args);
    
    const bool isHeterogeneous() const { return heterogeneous_;}
   
    void resetHeterogeneous() {  heterogeneous_=false; }
    
    void setHeterogeneous();
 
    const std::vector<bool> convergent() const { return convergent_; }

    const bool isDeputy();
  protected:
    std::vector<T>  args_;
    Ranges          ranges_;
    Params          params_;
    unsigned        distId_;
    bool heterogeneous_;
    std::vector<bool> convergent_;

  private:
    void extend (unsigned range_prod);

    void cartesianProduct (
      Params::const_iterator first2, Params::const_iterator last2);
    
    std::vector<bool> incrementLabel(std::vector<bool> label);
    
    bool lessLabels(std::vector<bool> labelA, std::vector<bool> labelB);
    
    bool orLabels(std::vector<bool> labelA, std::vector<bool> labelB, std::vector<bool> labelRes);
    
    std::vector<bool> resetLabel(std::vector<bool> label);
};

}  // namespace Horus

#endif  // YAP_PACKAGES_CLPBN_HORUS_GENERICFACTOR_H_

