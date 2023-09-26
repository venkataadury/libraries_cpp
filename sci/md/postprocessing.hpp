#include "coordinates.hpp"
#include <cassert>
#include <map>
namespace postprocessing
{
    using MD::numeric;
    using MD::Trajectory;
    using MD::Frame;

    class NoFramesInTrajectoryException : public std::exception {};

    struct DiffusionOutput
    {
      std::string component;
      numeric delt;
      numeric TMSD;
    };

    class DiffusionCoefficient
    {
        int skframe=0;
        numeric dt;
    public:
        DiffusionCoefficient(numeric t=1.0,int skip_frames=0)
        {
            this->dt=t;
            this->skframe=skip_frames;
        }

        std::vector<std::pair<std::string,numeric>> computeTMSD(const Trajectory& traj, int df, std::vector<int> atom_idxs=std::vector<int>(),int fstep=1) /*df is no. of frames gapped*/
        {
            assert(df>0);
            int nframes=traj.getNumFrames();
            if(skframe+df+1 >= nframes) throw NoFramesInTrajectoryException();
            std::map<std::string,numeric> records;
            bool setup=false;
            long Nfr=0;
            for(int i=skframe;i<traj.getNumFrames();i+=fstep)
            {
                const Frame& fr = traj.getFrame(i);
                if(!setup)
                {
                    setup=true;
                    if(!atom_idxs.size())
                    {
                        atom_idxs=std::vector<int>();
                        for(int k=0;k<fr.getNumAtoms();k++) atom_idxs.push_back(k);
                    }
                    for(int aidx : atom_idxs) records[fr.getAtom(aidx)->getName()]=0;
                }
                if(i+df>=nframes) break;
                const Frame& fr2=traj.getFrame(i+df);
                for(int aidx : atom_idxs) records[fr.getAtom(aidx)->getName()]+=fr.getAtom(aidx)->getDistance2To(fr2.getAtom(aidx));
                Nfr++;
            }
            std::vector<std::pair<std::string,numeric>> ret;
            for (const auto& [key, value] : records) ret.push_back(std::make_pair(key,value/Nfr));
            return ret;
        }

        std::vector<std::pair<std::string,numeric>> getDiffusionCoefficientAt(numeric dt,const Trajectory& traj,numeric timestep_per_frame,std::vector<int> atom_idxs=std::vector<int>(),int fstep=1)
        {
            int df=::round(dt/timestep_per_frame);
            std::vector<std::pair<std::string,numeric>> ret=computeTMSD(traj,df,atom_idxs,fstep);
            const numeric den=6*dt;
            std::map<std::string,int> Nat;

            std::vector<MD::Atom*> atoms=traj.getFrame(0).getAtoms();
            if(!atom_idxs.size())
            {
                atom_idxs=std::vector<int>();
                for(int k=0;k<atoms.size();k++) atom_idxs.push_back(k);
            }
            for(int idx : atom_idxs)
            {
                std::string nname=atoms[idx]->getName();
                if(Nat.find(nname)==Nat.end()) Nat[nname]=1;
                else Nat[nname]++;
            }
            for(auto& p : ret) p.second/=(Nat[p.first]*den); /*(traj-coord-units)^2/timestep-units*/
            return ret;
        }

        std::map<std::string,std::vector<numeric>> getDiffusionCoefficients(const Trajectory& traj,numeric timestep_per_frame,numeric dt_s,numeric dt_e,std::vector<int> atom_idxs=std::vector<int>(),int fstep=1)
        {
            int df_s = ::round(dt_s/timestep_per_frame);
            int df_e = ::round(dt_e/timestep_per_frame);
            assert(df_s<=df_e);
            numeric del_dt=timestep_per_frame; //*(df_e-df_s);
            std::cout <<"# Using "<<del_dt<<" as stride between different 'dt's for diffusion coefficient\n";

            std::map<std::string,std::vector<numeric>> ret;
            float dt=dt_s;
            while(dt<dt_e)
            {
                std::vector<std::pair<std::string,numeric>> cdc=getDiffusionCoefficientAt(dt,traj,timestep_per_frame,atom_idxs,fstep);
                for(const auto& p : cdc)
                {
                    if(ret.find(p.first)==ret.end()) ret[p.first]=std::vector<numeric>();
                    ret[p.first].push_back(p.second);
                }
                dt+=del_dt;
                //std::cout <<"# Reached dt="<<dt<<" of "<<dt_e<<"\n";
            }
            return ret;
        }

        void writeTable(const std::map<std::string,std::vector<numeric>>& tab,std::ostream& os,const std::string& index_name="SrNo",numeric step_index=1,numeric start_idx=0)
        {
            os << "# "<<index_name;
            for (const auto& [key, value] : tab) os << " "<<key;
            os<<"\n";
            numeric idx=start_idx;
            std::vector<std::vector<numeric>> matrix;
            for (const auto& [key, value] : tab) matrix.push_back(value);
            assert(matrix.size());
            for(int i=0;i<matrix[0].size();i++)
            {
                os << idx;
                for(int j=0;j<matrix.size();j++) os <<" "<< matrix[j][i];
                os << "\n";
                idx+=step_index;
            }
        }

        /*std::vector<DiffusionOutput> compute(Trajectory traj, int nstart=50, int nend=-1,int sstep=1,std::vector<int> atom_idxs=std::vector<int>(),int fstep=1)
        {
            std::vector<DiffusionOutput> ret;
            if(nend<0) nend=traj.getNumFrames()/2;
            for(int delt=nstart;delt<nend;delt+=sstep)
            {
                std::vector<std::pair<std::string,numeric>> tmsd_data=computeTMSD(traj,delt,atom_idxs,fstep);
                ret.push_back(DiffusionOutput({}));
            }
            return ret;
        }*/
    };
}
