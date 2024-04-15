#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
#if defined(__cplusplus)
extern "C" {
#endif

extern void _AMPA_syn_reg(void);
extern void _GABA_syn_reg(void);
extern void _kd_current_reg(void);
extern void _leak_current_reg(void);
extern void _na_current_reg(void);

void modl_reg() {
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");
    fprintf(stderr, " \"AMPA_syn.mod\"");
    fprintf(stderr, " \"GABA_syn.mod\"");
    fprintf(stderr, " \"kd_current.mod\"");
    fprintf(stderr, " \"leak_current.mod\"");
    fprintf(stderr, " \"na_current.mod\"");
    fprintf(stderr, "\n");
  }
  _AMPA_syn_reg();
  _GABA_syn_reg();
  _kd_current_reg();
  _leak_current_reg();
  _na_current_reg();
}

#if defined(__cplusplus)
}
#endif
