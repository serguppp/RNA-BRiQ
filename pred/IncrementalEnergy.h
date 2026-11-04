#ifndef BRIQ_INCREMENTAL_ENERGY_H
#define BRIQ_INCREMENTAL_ENERGY_H

#include "pred/BRFoldingTree.h"
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace NSPpred {

// Disjoint sums matching the PDB-ownership pairs.csv/nonpair.csv exporter.
struct EnergySums {
    double bb=0, bo=0, oo=0, clash=0;
    double baseRibose=0, basePhosphate=0, ribosePhosphate=0, rotamer=0, phosphate=0;
    double pairTotal() const { return bb+bo+oo+clash; }
    double nonpairTotal() const { return baseRibose+basePhosphate+ribosePhosphate+rotamer+phosphate; }
    double total() const { return pairTotal()+nonpairTotal(); }
    std::array<double,9> values() const {
        return {{bb,bo,oo,clash,baseRibose,basePhosphate,ribosePhosphate,rotamer,phosphate}};
    }
    bool agrees(double energy) const {
        for (double value : values()) if (!std::isfinite(value)) return false;
        return std::isfinite(energy) && std::fabs(total()-energy)<=1e-9+1e-12*std::fabs(energy);
    }
    void interaction(int first, int second, double oxygen, double steric, int kind) {
        if (first<0 || second<0) {
            if (oxygen!=0 || steric!=0) throw std::logic_error("Energy assigned to absent phosphate");
        } else if (first==second) {
            if (kind==0) baseRibose+=oxygen+steric;
            else if (kind==1) basePhosphate+=oxygen+steric;
            else if (kind==2) ribosePhosphate+=oxygen+steric;
            else throw std::logic_error("Unexpected local phosphate-phosphate energy");
        } else {
            if (kind<2) bo+=oxygen; else oo+=oxygen;
            clash+=steric;
        }
    }
};

// Cache only numerically identical node states. No geometric tolerance or
// distance approximation is introduced; every changed node invalidates all
// of its pairwise and directed interactions, including phosphate neighbors.
class IncrementalEnergy {
    std::vector<std::vector<double> > previous;
    std::vector<int> previousSeparation;
    std::vector<std::array<double,4> > pairs;
    std::vector<std::array<double,3> > directed;
    std::vector<std::array<double,2> > riboseParts;
    std::vector<std::array<double,6> > directedParts;
    bool detailsReady=false;
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

    double score(BRFoldingTree& tree, EnergySums* sums=NULL) {
        const size_t n=tree.seqLen;
        const bool reset=previous.size()!=n || (sums!=NULL && !detailsReady);
        if (reset) {
            pairs.resize(n*n);
            directed.resize(n*n);
            if (sums!=NULL) {riboseParts.resize(n*n);directedParts.resize(n*n);}
        }
        if (sums!=NULL) *sums=EnergySums();
        auto phosphateOwner=[&](size_t i) {
            return tree.nodes[i]->connectToNeighbor && i+1<n ? static_cast<int>(i+1) : -1;
        };
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
                double oxygen=0,steric=0;
                pairs[index]={{getBaseBaseEnergy(a,b,sep,tree.et,false),
                    baseBaseClash(a,b,sep,tree.et,0.0,false),
                    getRiboseRiboseEnergy(a,b,sep,tree.et,false,sums ? &oxygen : NULL,sums ? &steric : NULL),
                    getPhoPhoEnergy(a,b,sep,tree.et,false)}};
                if (sums!=NULL) riboseParts[index]={{oxygen,steric}};
            }
            // Preserve the full scorer's accumulation order and total exactly.
            for (size_t k=0;k<4;++k) total+=pairs[index][k];
            if (sums!=NULL) {
                sums->bb+=pairs[index][0];
                sums->oo+=riboseParts[index][0];
                sums->clash+=pairs[index][1]+riboseParts[index][1];
                // Two distinct connected nodes have distinct phosphate owners.
                sums->interaction(phosphateOwner(i),phosphateOwner(j),0,pairs[index][3],3);
            }
        }
        for (size_t i=0;i<n;++i) for (size_t j=0;j<n;++j) {
            const size_t index=i*n+j;
            const int sep=tree.sepTable[index];
            if (reset || changed[i] || changed[j] || sep!=previousSeparation[index]) {
                ++directedEvaluations;
                BRNode *a=tree.nodes[i], *b=tree.nodes[j];
                std::array<double,6> detail={{0,0,0,0,0,0}};
                directed[index]={{getBaseRiboseEnergy(a,b,sep,tree.et,false,sums ? &detail[0] : NULL,sums ? &detail[1] : NULL),
                    getBasePhoEnergy(a,b,sep,tree.et,false,sums ? &detail[2] : NULL,sums ? &detail[3] : NULL),
                    getRibosePhoEnergy(a,b,sep,tree.et,false,sums ? &detail[4] : NULL,sums ? &detail[5] : NULL)}};
                if (sums!=NULL) directedParts[index]=detail;
            }
            for (size_t k=0;k<3;++k) total+=directed[index][k];
            if (sums!=NULL) {
                const auto& d=directedParts[index];
                sums->interaction(i,j,d[0],d[1],0);
                sums->interaction(i,phosphateOwner(j),d[2],d[3],1);
                sums->interaction(i,phosphateOwner(j),d[4],d[5],2);
            }
        }
        for (size_t i=0;i<n;++i) {
            total+=tree.nodes[i]->rot->energy;
            if (sums!=NULL) sums->rotamer+=tree.nodes[i]->rot->energy;
        }
        for (size_t i=0;i<n;++i) {
            total+=tree.nodes[i]->pho.e;
            if (sums!=NULL) sums->phosphate+=tree.nodes[i]->pho.e;
        }
        previous.swap(current);
        previousSeparation.assign(tree.sepTable,tree.sepTable+n*n);
        detailsReady=sums!=NULL;
        if (sums!=NULL && !sums->agrees(total)) throw std::runtime_error("Energy component sums do not reproduce total");
        return total;
    }
};
}
#endif
