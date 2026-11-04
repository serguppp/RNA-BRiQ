#ifndef BRIQ_INCREMENTAL_ENERGY_H
#define BRIQ_INCREMENTAL_ENERGY_H

#include "pred/BRFoldingTree.h"
#include <array>
#include <vector>

namespace NSPpred {

// Cache only numerically identical node states. No geometric tolerance or
// distance approximation is introduced; every changed node invalidates all
// of its pairwise and directed interactions, including phosphate neighbors.
class IncrementalEnergy {
    std::vector<std::vector<double> > previous;
    std::vector<int> previousSeparation;
    std::vector<std::array<double,4> > pairs;
    std::vector<std::array<double,3> > directed;
    static void append(std::vector<double>& out, const XYZ& p) {
        out.push_back(p.x_); out.push_back(p.y_); out.push_back(p.z_);
    }
    static std::vector<double> state(BRNode* node) {
        std::vector<double> out;
        out.reserve(130);
        out.push_back(node->seqID);
        out.push_back(node->baseType);
        out.push_back(node->baseAtomNum);
        out.push_back(node->connectToNeighbor);
        append(out, node->cs1.origin_);
        for (int i=0;i<3;++i) for (int j=0;j<3;++j) {
            out.push_back(node->cs1.localToGlobalTM.mtx[i][j]);
            out.push_back(node->cs1.globalToLocalTM.mtx[i][j]);
        }
        for (int i=0;i<node->baseAtomNum;++i) append(out,node->baseAtomCoords[i]);
        for (int i=0;i<11;++i) append(out,node->riboAtomCoords[i]);
        for (int i=0;i<7;++i) append(out,node->pho.tList[i]);
        out.push_back(node->rot->energy);
        out.push_back(node->pho.e);
        return out;
    }
 public:
    size_t changedNodes = 0;
    size_t pairEvaluations = 0;
    size_t directedEvaluations = 0;

    double score(BRFoldingTree& tree) {
        const size_t n=tree.seqLen;
        const bool reset=previous.size()!=n;
        if (reset) {
            pairs.resize(n*n);
            directed.resize(n*n);
        }
        std::vector<std::vector<double> > current;
        std::vector<bool> changed(n);
        changedNodes=pairEvaluations=directedEvaluations=0;
        for (size_t i=0;i<n;++i) {
            current.push_back(state(tree.nodes[i]));
            changed[i]=reset || current[i]!=previous[i];
            changedNodes+=changed[i];
        }
        double total=0;
        for (size_t i=0;i<n;++i) for (size_t j=i+1;j<n;++j) {
            const size_t index=i*n+j;
            const int sep=tree.sepTable[index];
            if (reset || changed[i] || changed[j] || sep!=previousSeparation[index]) {
                ++pairEvaluations;
                BRNode *a=tree.nodes[i], *b=tree.nodes[j];
                pairs[index]={{getBaseBaseEnergy(a,b,sep,tree.et,false),
                    baseBaseClash(a,b,sep,tree.et,0.0,false),
                    getRiboseRiboseEnergy(a,b,sep,tree.et,false),
                    getPhoPhoEnergy(a,b,sep,tree.et,false)}};
            }
            // Preserve the full scorer's accumulation order.
            for (size_t k=0;k<4;++k) total+=pairs[index][k];
        }
        for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) {
            const size_t index=i*n+j;
            const int sep=tree.sepTable[index];
            if (reset || changed[i] || changed[j] || sep!=previousSeparation[index]) {
                ++directedEvaluations;
                BRNode *a=tree.nodes[i], *b=tree.nodes[j];
                directed[index]={{getBaseRiboseEnergy(a,b,sep,tree.et,false),
                    getBasePhoEnergy(a,b,sep,tree.et,false),
                    getRibosePhoEnergy(a,b,sep,tree.et,false)}};
            }
            for (size_t k=0;k<3;++k) total+=directed[index][k];
        }
        for (size_t i=0;i<n;++i) total+=tree.nodes[i]->rot->energy;
        for (size_t i=0;i<n;++i) total+=tree.nodes[i]->pho.e;
        previous.swap(current);
        previousSeparation.assign(tree.sepTable,tree.sepTable+n*n);
        return total;
    }
};
}
#endif
