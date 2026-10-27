// Test-only oracle. Compile against the unmodified tzok commit, never the exporter.
#include "pred/BRFoldingTree.h"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cstdlib>
using namespace NSPpred;
using namespace std;

int main(int argc, char** argv) {
    if(argc != 3 && argc != 6) return 2;
    srand(0);
    EnergyTable* et = new EnergyTable();
    BRFoldingTree* ft = new BRFoldingTree(argv[1], et);
    double bw=1, cw=1, aw=1;
    if(argc == 6) { bw=atof(argv[3]); cw=atof(argv[4]); aw=atof(argv[5]); }
    ofstream out(argv[2]);
    out << setprecision(17);
    out << "kind\ti\tj\tsep\tenergy\toxygen\tclash\n";
    auto emit = [&](const char* kind,int i,int j,int sep,double e,double oxygen,double clash) {
        out << kind << '\t' << i << '\t' << j << '\t' << sep << '\t' << e << '\t' << oxygen << '\t' << clash << '\n';
    };
    // The helper return values below are from the original, unchanged header.
    for(int i=0;i<ft->seqLen;i++) for(int j=i+1;j<ft->seqLen;j++) {
        BRNode* a=ft->nodes[i]; BRNode* b=ft->nodes[j]; int sep=ft->sepTable[i*ft->seqLen+j];
        double bb=getBaseBaseEnergy(a,b,sep,et,false);
        double bc=baseBaseClash(a,b,sep,et,0,false)*aw;
        double rr=getRiboseRiboseEnergy(a,b,sep,et,false)*aw;
        double rh=0;
        if(squareDistance(a->riboAtomCoords[2],b->riboAtomCoords[2])<81 && abs(sep)>=2)
            rh=et->atET->getRiboseRiboseHbondEnergy(a->riboAtomCoords,b->riboAtomCoords)*aw;
        double pp=getPhoPhoEnergy(a,b,sep,et,false)*aw;
        emit("base_base",i,j,sep,bb,0,0);
        emit("base_clash",i,j,sep,bc,0,bc);
        emit("ribose_clash",i,j,sep,rr,rh,rr-rh);
        emit("phosphate_clash",i,j,sep,pp,0,pp);
    }
    for(int i=0;i<ft->seqLen;i++) for(int j=0;j<ft->seqLen;j++) {
        BRNode* a=ft->nodes[i]; BRNode* b=ft->nodes[j]; int sep=ft->sepTable[i*ft->seqLen+j];
        double br=getBaseRiboseEnergy(a,b,sep,et,false)*aw;
        double bp=getBasePhoEnergy(a,b,sep,et,false)*aw;
        double rp=getRibosePhoEnergy(a,b,sep,et,false)*aw;
        double boR=0,boP=0,roP=0;
        if(sep!=0 && squareDistance(a->baseAtomCoords[0],b->riboAtomCoords[2])<144) {
            const int core[]={5,6,4}, support[]={8,9,10};
            for(int k=0;k<3;k++) boR+=et->roET->getEnergy(a->baseType,k,global2local(a->cs1,b->riboAtomCoords[core[k]]),global2local(a->cs1,b->riboAtomCoords[support[k]]),sep);
        }
        if(sep!=0 && b->connectToNeighbor) {
            for(int k=1;k<=3;k++) boP+=et->roET->getEnergy(a->baseType,k==1?3:4,global2local(a->cs1,b->pho.tList[k]),global2local(a->cs1,b->pho.tList[k+3]),sep);
        }
        if(sep>2 && b->connectToNeighbor) roP=et->atET->getRibosePhoHbondEnergy(a->riboAtomCoords,b->pho.tList);
        boR*=et->para.wtBaseOxygen*aw; boP*=et->para.wtBaseOxygen*aw; roP*=aw;
        emit("base_ribose",i,j,sep,br,boR,br-boR);
        emit("base_phosphate",i,j,sep,bp,boP,bp-boP);
        emit("ribose_phosphate",i,j,sep,rp,roP,rp-roP);
    }
    for(int i=0;i<ft->seqLen;i++) {
        emit("rotamer",i,i,0,ft->nodes[i]->rot->energy,0,0);
        emit("phosphate_local",i,i,0,ft->nodes[i]->pho.e*(ft->chainBreakPoints[i]?bw:cw),0,0);
        // Node i owns the outgoing phosphate link; PDB atom ownership is i+1.
        emit("phosphate_owner",i,ft->connectToDownstream[i]?i+1:-1,0,0,0,0);
    }
    emit("global",-1,-1,0,ft->totalEnergy(bw,cw,aw,0,false),0,0);
    ofstream oxygen((string(argv[2])+".oxygen.tsv").c_str());
    oxygen << setprecision(17) << "i\tj\ttype\tsep\tvalue\treason\n";
    auto inspectOxygen = [&](int i,int j,int type,const XYZ& core,const XYZ& support,int sep) {
        BRNode* a=ft->nodes[i]; XYZ p=global2local(a->cs1,core);
        double value=et->roET->getEnergy(a->baseType,type,p,global2local(a->cs1,support),sep);
        string reason="nonzero";
        if(value==0) {
            if(sep==0 || (sep==1 && type!=4) || (sep==-1 && type!=2 && type!=4)) reason="separation_exclusion";
            else if(p.x_<=-3 || p.x_>=12 || p.y_<=-6 || p.y_>=10 || p.z_<=-5 || p.z_>=5) reason="outside_grid";
            else {
                int index=int((p.x_+3)*5)*4000+int((p.y_+6)*5)*50+int((p.z_+5)*5);
                double raw=sep==1?et->roET->etListP1[a->baseType][index]:sep==-1?et->roET->etListM1[a->baseType*2+(type==4)][index]:et->roET->etList[a->baseType*5+type][index];
                reason=raw==0?"zero_table_entry":"angular_rejection";
            }
        }
        oxygen << i << '\t' << j << '\t' << type << '\t' << sep << '\t' << value << '\t' << reason << '\n';
    };
    for(int i=0;i<ft->seqLen;i++) for(int j=0;j<ft->seqLen;j++) {
        BRNode* a=ft->nodes[i]; BRNode* b=ft->nodes[j]; int sep=ft->sepTable[i*ft->seqLen+j];
        if(sep!=0 && squareDistance(a->baseAtomCoords[0],b->riboAtomCoords[2])<144) {
            const int core[]={5,6,4}, support[]={8,9,10};
            for(int k=0;k<3;k++) inspectOxygen(i,j,k,b->riboAtomCoords[core[k]],b->riboAtomCoords[support[k]],sep);
        }
        if(sep!=0 && b->connectToNeighbor)
            for(int k=1;k<=3;k++) inspectOxygen(i,j,k==1?3:4,b->pho.tList[k],b->pho.tList[k+3],sep);
    }
    ofstream hb((string(argv[2])+".hbond.tsv").c_str());
    hb << setprecision(17) << "kind\ti\tj\tsep\tdistance\tvalue\treason\n";
    for(int i=0;i<ft->seqLen;i++) for(int j=0;j<ft->seqLen;j++) {
        BRNode* a=ft->nodes[i]; BRNode* b=ft->nodes[j]; int sep=ft->sepTable[i*ft->seqLen+j];
        if(i<j) {
            double distance=a->riboAtomCoords[5].distance(b->riboAtomCoords[5]);
            bool allowed=abs(sep)>=2 && squareDistance(a->riboAtomCoords[2],b->riboAtomCoords[2])<81;
            double value=allowed?et->atET->getRiboseRiboseHbondEnergy(a->riboAtomCoords,b->riboAtomCoords):0;
            hb << "O2O2\t" << i << '\t' << j << '\t' << sep << '\t' << distance << '\t' << value << '\t'
               << (!allowed?"separation_or_outer_cutoff":distance>=4?"distance_cutoff":value==0?"zero_table_entry":"nonzero") << '\n';
        }
        if(b->connectToNeighbor) {
            double distance=min(a->riboAtomCoords[5].distance(b->pho.tList[2]),a->riboAtomCoords[5].distance(b->pho.tList[3]));
            bool allowed=sep>2;
            double value=allowed?et->atET->getRibosePhoHbondEnergy(a->riboAtomCoords,b->pho.tList):0;
            hb << "O2OP\t" << i << '\t' << j << '\t' << sep << '\t' << distance << '\t' << value << '\t'
               << (!allowed?"separation_exclusion":distance>=4?"distance_cutoff":value==0?"zero_table_entry":"nonzero") << '\n';
        }
    }
    ofstream parameters((string(argv[2])+".parameters.json").c_str());
    parameters << setprecision(17) << "{\"base_oxygen_weight\":" << et->para.wtBaseOxygen << ",\"clash_weight\":" << aw << "}";
    // One-body/local terms are anchored at a node, not necessarily single-residue physics.
    cerr << setprecision(17) << "global " << ft->totalEnergy(bw,cw,aw,0,false) << '\n';
    // Original ownership/destructor behavior is retained in the reference build.
    delete ft; delete et;
    return out.good() && oxygen.good()?0:1;
}
