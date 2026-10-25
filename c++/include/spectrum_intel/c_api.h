/* spectrum_intel C API: for embedding the C++ library from C, Go, Rust, Python (ctypes)...
 *
 *   si_engine* e = si_engine_create("{\"disable\": [\"dji\"]}", &err);
 *   char* json = si_analyze(e, samples, n, SI_CI16, 61.44e6, 2.437e9, 56e6,
 *                           "{\"start_time_utc\": \"2026-10-09T10:00:00Z\"}");
 *   ...  json is the response (schema/response.schema.json)
 *   si_free(json);
 *   si_engine_destroy(e);
 */
#ifndef SPECTRUM_INTEL_C_API_H
#define SPECTRUM_INTEL_C_API_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct si_engine si_engine;

enum si_iq_format { SI_CF32 = 0, SI_CF64 = 1, SI_CI16 = 2, SI_CI8 = 3, SI_CU8 = 4 };

/* config_json: NULL or {"enable": [...], "disable": [...], "threads": N, "compact": bool}.
 * Returns NULL on error; *error (if not NULL) then receives a message to si_free(). */
si_engine* si_engine_create(const char* config_json, char** error);
void si_engine_destroy(si_engine* engine);

/* bandwidth_hz <= 0: the whole sampled band. meta_json: NULL or {"start_time_utc", "receiver_ppm",
 * "gain_db", "source_id", "request_id"}. Returns the response JSON (never NULL); free with si_free(). */
char* si_analyze(const si_engine* engine, const void* iq, size_t n_samples, int format, double sample_rate_hz,
                 double center_frequency_hz, double bandwidth_hz, const char* meta_json);

/* schema/requirements.json as compiled in; free with si_free(). */
char* si_requirements(void);
const char* si_version(void);
void si_free(char* p);

#ifdef __cplusplus
}
#endif

#endif /* SPECTRUM_INTEL_C_API_H */
