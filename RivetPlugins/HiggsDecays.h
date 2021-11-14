#ifndef _HIGGSDECAYS_H_
#define _HIGGSDECAYS_H_

#include "Rivet/Analysis.hh"

#include "Rivet/Tools/RivetYODA.hh"

namespace Rivet
{
    class HiggsDecays : public Analysis 
    {
    public:
    //---ctors---
        HiggsDecays():
            Analysis("HiggsDecays"),
            sumW_(0.)
            {};

    //---dtor---
        ~HiggsDecays() {};

    //---utils--
        void init();
        void analyze(const Event& event);
        void finalize();

    private:
        double     sumW_;
        Histo1DPtr h_sumW_;
    };

    DECLARE_RIVET_PLUGIN(HiggsDecays);
}

#endif
