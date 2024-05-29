
#ifndef NEURODEGENERATION_MODEL_HPP
#define NEURODEGENERATION_MODEL_HPP

#include <deal.II/base/point.h>

using namespace dealii;

class NeurodegenerationModel {
public:
    enum class DiseaseType { AmyloidBeta, Tau, AlphaSynuclein, TDP43 };
    
    // Costruttore predefinito
    NeurodegenerationModel() : disease(DiseaseType::AmyloidBeta) {}

    // Nuovo costruttore che accetta un DiseaseType
    NeurodegenerationModel(DiseaseType disease_) : disease(disease_) {}

    //void set_disease(DiseaseType disease_) { disease = disease_; }
    
    bool is_in_disease_region(const Point<3> &p) const;
    bool is_in_disease_region(const Point<2> &p) const;
private:
    DiseaseType disease;
    
};


#endif // NEURODEGENERATION_MODEL_HPP
