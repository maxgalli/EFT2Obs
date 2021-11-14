#include "HiggsDecays.h"

namespace Rivet
{
    void HiggsDecays::init()
    {
        //---Histograms
        book(h_sumW_, "sumWeights", {-1, 1});
    }

    void HiggsDecays::analyze(const Event& event)
    {               
        h_sumW_->fill(0, event.weights()[0]);
        return;
    }

    void HiggsDecays::finalize()
    {
        return;
    }
}
