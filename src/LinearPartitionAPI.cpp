#include "LinearPartitionAPI.h"
#include "LinearPartition.h"

#include <algorithm>
#include <cctype>
#include <exception>
#include <mutex>
#include <string>
#include <cstring>

namespace {

thread_local std::string g_last_error;

int sanitize_parameters(int beamsize) {
    return beamsize > 0 ? beamsize : 100;
}

std::string sanitize_sequence(const char* sequence) {
    std::string seq(sequence);
    if (seq.empty()) {
        throw std::invalid_argument("Sequence is empty");
    }
    std::transform(seq.begin(), seq.end(), seq.begin(), [](unsigned char c) {
        char upper = static_cast<char>(std::toupper(c));
        return upper == 'T' ? 'U' : upper;
    });
    return seq;
}

void set_error(const std::string& message) {
    g_last_error = message;
}

}  // namespace

extern "C" {

const char* lp_get_last_error(void) {
    return g_last_error.c_str();
}

int lp_compute_ensemble_energy(const char* sequence,
                               int beamsize,
                               int no_sharp_turn,
                               int verbose,
                               int dangles,
                               double* out_energy) {
    if (sequence == nullptr || out_energy == nullptr) {
        set_error("Null pointer argument");
        return LP_STATUS_INVALID_ARGUMENT;
    }

    try {
        std::string seq = sanitize_sequence(sequence);
        const int sanitized_beam = sanitize_parameters(beamsize);

        BeamCKYParser parser(
            sanitized_beam,
            no_sharp_turn != 0,
            verbose != 0,
            /*bppfile=*/"",
            /*bppfileindex=*/"",
            /*pf_only=*/true,
            /*bpp_cutoff=*/0.0f,
            /*forestfile=*/"",
            /*mea_=*/false,
            /*gamma=*/3.0f,
            /*mea_file_index=*/"",
            /*bpseq=*/false,
            /*threshknot_=*/false,
            /*threshknot_threshold=*/0.3f,
            /*threshknot_file_index=*/"",
            /*shape_file_path=*/"",
            /*is_fasta=*/false,
            /*dangles=*/dangles > 0 ? dangles : 2,
            /*quiet=*/true);

        double energy = parser.parse(seq);
        *out_energy = energy;
        set_error("");
        return LP_STATUS_SUCCESS;
    } catch (const std::exception& ex) {
        set_error(ex.what());
        return LP_STATUS_EXCEPTION;
    } catch (...) {
        set_error("Unknown error");
        return LP_STATUS_EXCEPTION;
    }
}

int lp_compute_bpp_matrix(const char* sequence,
                          int beamsize,
                          int no_sharp_turn,
                          int verbose,
                          int dangles,
                          double cutoff,
                          double* out_matrix,
                          long long out_matrix_len) {
    if (sequence == nullptr || out_matrix == nullptr) {
        set_error("Null pointer argument");
        return LP_STATUS_INVALID_ARGUMENT;
    }

    try {
        std::string seq = sanitize_sequence(sequence);
        const int sanitized_beam = sanitize_parameters(beamsize);
        const size_t n = seq.size();
        const size_t required = n * n;
        if (out_matrix_len < static_cast<long long>(required)) {
            set_error("Output matrix buffer too small");
            return LP_STATUS_INVALID_ARGUMENT;
        }

        BeamCKYParser parser(
            sanitized_beam,
            no_sharp_turn != 0,
            verbose != 0,
            /*bppfile=*/"",
            /*bppfileindex=*/"",
            /*pf_only=*/false,
            /*bpp_cutoff=*/static_cast<float>(cutoff),
            /*forestfile=*/"",
            /*mea_=*/false,
            /*gamma=*/3.0f,
            /*mea_file_index=*/"",
            /*bpseq=*/false,
            /*threshknot_=*/false,
            /*threshknot_threshold=*/0.3f,
            /*threshknot_file_index=*/"",
            /*shape_file_path=*/"",
            /*is_fasta=*/false,
            /*dangles=*/dangles > 0 ? dangles : 2,
            /*quiet=*/true);

        (void)parser.parse(seq);
        parser.export_bpp_dense(out_matrix, static_cast<size_t>(out_matrix_len));
        set_error("");
        return LP_STATUS_SUCCESS;
    } catch (const std::exception& ex) {
        set_error(ex.what());
        return LP_STATUS_EXCEPTION;
    } catch (...) {
        set_error("Unknown error");
        return LP_STATUS_EXCEPTION;
    }
}

int lp_calculate_mea_structure(const char* sequence,
                               int beamsize,
                               int no_sharp_turn,
                               int verbose,
                               int dangles,
                               double gamma,
                               char* out_structure,
                               long long out_structure_len,
                               double* out_energy) {
    if (sequence == nullptr || out_structure == nullptr || out_energy == nullptr) {
        set_error("Null pointer argument");
        return LP_STATUS_INVALID_ARGUMENT;
    }

    try {
        std::string seq = sanitize_sequence(sequence);
        const int sanitized_beam = sanitize_parameters(beamsize);
        const size_t n = seq.size();
        if (out_structure_len < static_cast<long long>(n + 1)) {
            set_error("Structure buffer too small");
            return LP_STATUS_INVALID_ARGUMENT;
        }

        BeamCKYParser parser(
            sanitized_beam,
            no_sharp_turn != 0,
            verbose != 0,
            /*bppfile=*/"",
            /*bppfileindex=*/"",
            /*pf_only=*/false,
            /*bpp_cutoff=*/0.0f,
            /*forestfile=*/"",
            /*mea_=*/true,
            /*gamma=*/static_cast<float>(gamma),
            /*mea_file_index=*/"",
            /*bpseq=*/false,
            /*threshknot_=*/false,
            /*threshknot_threshold=*/0.3f,
            /*threshknot_file_index=*/"",
            /*shape_file_path=*/"",
            /*is_fasta=*/false,
            /*dangles=*/dangles > 0 ? dangles : 2,
            /*quiet=*/true);

        double energy = parser.parse(seq);
        const std::string& structure = parser.get_last_mea_structure();
        if (structure.empty()) {
            set_error("MEA structure not available");
            return LP_STATUS_EXCEPTION;
        }
        std::strncpy(out_structure, structure.c_str(), static_cast<size_t>(out_structure_len));
        out_structure[static_cast<size_t>(n)] = '\0';
        *out_energy = energy;
        set_error("");
        return LP_STATUS_SUCCESS;
    } catch (const std::exception& ex) {
        set_error(ex.what());
        return LP_STATUS_EXCEPTION;
    } catch (...) {
        set_error("Unknown error");
        return LP_STATUS_EXCEPTION;
    }
}

}  // extern "C"

