#include "neurodegeneration_models.hpp"

bool NeurodegenerationModel::is_in_disease_region(const Point<3> &p) const {
    switch (disease) {
        case DiseaseType::AmyloidBeta:
            return (p[0] >= 45.0 && p[0] <= 83.0 && p[1] >= 85.0 && p[1] <= 110.0 && p[2] >= 90.0 && p[2] <= 116.0) ||
                   (p[0] >= 45.0 && p[0] <= 83.0 && p[1] >= 40.0 && p[1] <= 150.0 && p[2] >= 11.0 && p[2] <= 50.0);
        case DiseaseType::Tau:
            return p[0] >= 80.0 && p[0] <= 83.0 && p[1] >= 53.0 && p[1] <= 55.0 && p[2] >= 63.0 && p[2] <= 65.0;
        case DiseaseType::AlphaSynuclein:
            return p[0] >= 45.0 && p[0] <= 63.0 && p[1] >= 40.0 && p[1] <= 60.0 && p[2] >= 11.0 && p[2] <= 35.0;
        case DiseaseType::TDP43:
            return (p[0] >= 45.0 && p[0] <= 63.0 && p[1] >= 45.0 && p[1] <= 55.0 && p[2] >= 11.0 && p[2] <= 15.0) ||
                   (p[0] >= 45.0 && p[0] <= 83.0 && p[1] >= 65.0 && p[1] <= 85.0 && p[2] >= 80.0 && p[2] <= 116.0);
        default:
            return false;
    }
}

// Implementazione per i punti a 2 dimensioni
bool NeurodegenerationModel::is_in_disease_region(const Point<2> &p) const {
    // For 2D, we can ignore the z-coordinate checks or provide simplified regions.
        switch (disease) {
            case DiseaseType::AmyloidBeta:
                return (p[0] >= 45.0 && p[0] <= 83.0 && p[1] >= 85.0 && p[1] <= 110.0) ||
                       (p[0] >= 45.0 && p[0] <= 83.0 && p[1] >= 40.0 && p[1] <= 150.0);
            case DiseaseType::Tau:
                return p[0] >= 80.0 && p[0] <= 83.0 && p[1] >= 53.0 && p[1] <= 55.0;
            case DiseaseType::AlphaSynuclein:
                return p[0] >= 45.0 && p[0] <= 63.0 && p[1] >= 40.0 && p[1] <= 60.0;
            case DiseaseType::TDP43:
                return (p[0] >= 45.0 && p[0] <= 63.0 && p[1] >= 45.0 && p[1] <= 55.0) ||
                       (p[0] >= 45.0 && p[0] <= 83.0 && p[1] >= 65.0 && p[1] <= 85.0);
            default:
                return false;
        }
}