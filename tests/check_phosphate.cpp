// Exercise phosphate rebuilding without the multi-gigabyte base-pair tables.
#include "forcefield/PO3Builder.h"
#include "model/BaseRotamer.h"
#include "model/ProteinRep.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace NSPforcefield;
int main(int argc, char** argv) {
    if (argc != 2) return 2;
    Parameter parameters;
    PO3Builder builder(&parameters, NULL, NULL);
    RNAPDB pdb(argv[1], "test");
    const auto bases=pdb.getBaseList();
    size_t checked=0;
    for (size_t i=0; i+1<bases.size(); ++i) {
        RNABase *a=bases[i], *b=bases[i+1];
        if (!a->backboneComplete() || !b->backboneComplete()) continue;
        BaseRotamer ra(a), rb(b);
        LocalFrame cs=a->getCoordSystem()+ra.mv12;
        XYZ ribose[11];
        ribose[3]=b->getAtom("C4'")->coord;
        ribose[4]=b->getAtom("O4'")->coord;
        ribose[7]=b->getAtom("C5'")->coord;
        for (int sign=-1; sign<=1; sign+=2) {
            auto result=builder.getPhoLocal(cs,ribose,sign*std::abs(ra.improper),rb.improper);
            if (!std::isfinite(result.e)) throw std::runtime_error("Nonfinite phosphate");
            ++checked;
        }
    }
    if (!checked) return 1;
    std::cout << "Phosphate probes PASS: " << checked << std::endl;
}
