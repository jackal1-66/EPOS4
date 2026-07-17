//
//  This file is part of EPOS4
//  Copyright (C) 2022 research institutions and authors (See CREDITS file)
//  This file is distributed under the terms of the GNU General Public License version 3 or later
//  (See COPYING file for the text of the licence)
//

#include <math.h>
#include <iomanip>
#include <ctime>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <fstream>
#include <charconv>
#include <string>
#include "es.h"
#include "eo3.h"

using namespace std ;

extern void mix(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p);

double sign(double x)
{
        if(x>0) return 1. ;
        else if(x<0.) return -1. ;
        else return 0. ;
}


EoS3f::EoS3f(char *filename, double _B, double _volex0, double _delta0, double _aaa, double _bbb)
{
        // Bulk-read the table and tokenize with std::from_chars instead of
        // per-token ifstream>>. The istream/locale overhead is removed,
        // meaning tens of seconds per job on the ~440 MB eos4f table.
        string buf ;
        {
                ifstream fraw (filename, ios::binary) ;
                fraw.seekg(0, ios::end) ;
                streamoff len = fraw.tellg() ;
                if(len < 0) len = 0 ;
                buf.resize((size_t)len) ;
                fraw.seekg(0, ios::beg) ;
                if(len > 0) fraw.read(&buf[0], len) ;
        }
        // Binary table produced at build time by src/MS/eosconvert.cpp from ASCII
        // layout: magic, int ne,nn, double B,volex0,delta0,
        // aaa,bbb,emax,e0, egrid[ne], ngrid[ne*nn], pre/T/mub/muq/mus[npts].
        // Falls through to the ASCII parser when the `magic` is absent.
        static const char eosb_magic[16] = "EPOS4EOS3F-BIN1" ;
        if(buf.size() >= 16+2*sizeof(int)+7*sizeof(double)
           && memcmp(buf.data(), eosb_magic, 16) == 0){
                const char *bp = buf.data() + 16 ;
                int dims[2] ;
                double hdr[7] ;
                memcpy(dims, bp, sizeof dims) ; bp += sizeof dims ;
                memcpy(hdr, bp, sizeof hdr) ; bp += sizeof hdr ;
                ne = dims[0] ; nn = dims[1] ;
                B = hdr[0] ; volex0 = hdr[1] ; delta0 = hdr[2] ; aaa = hdr[3] ; bbb = hdr[4] ;
                emax = hdr[5] ; e0 = hdr[6] ;
        if( B > 0 )
        {
        if(fabs(B-_B)          >1e-5){ cout <<"\n\n\n      B =      "<<B     <<" instead of "<<_B     <<"; exiting\n\n\n\n" ; exit(0) ; }
        if(fabs(volex0-_volex0)>1e-5){ cout <<"\n\n\n      volex0 = "<<volex0<<" instead of "<<_volex0<<"; exiting\n\n\n\n" ; exit(0) ; }
        if(fabs(delta0-_delta0)>1e-5){ cout <<"\n\n\n      delta0 = "<<delta0<<" instead of "<<_delta0<<"; exiting\n\n\n\n" ; exit(0) ; }
        if(fabs(aaa-_aaa)      >1e-5){ cout <<"\n\n\n      aaa =    "<<aaa   <<" instead of "<<_aaa   <<"; exiting\n\n\n\n" ; exit(0) ; }
        if(fabs(bbb-_bbb)      >1e-5){ cout <<"\n\n\n      bbb =    "<<bbb   <<" instead of "<<_bbb   <<"; exiting\n\n\n\n" ; exit(0) ; }
        }
                const size_t npts = (size_t)ne*nn*nn*nn ;
                const size_t need = 16 + 2*sizeof(int)
                  + (7 + (size_t)ne + (size_t)ne*nn + 5*npts)*sizeof(double) ;
                if(ne <= 0 || nn <= 0 || buf.size() != need){
                        cout<<"reading failed, truncated or corrupt binary eos table"<<endl;
                        exit(13);
                }
                egrid = new double [ne] ;
                ngrid = new double [ne*nn] ;
                memcpy(egrid, bp, (size_t)ne*sizeof(double)) ;       bp += (size_t)ne*sizeof(double) ;
                memcpy(ngrid, bp, (size_t)ne*nn*sizeof(double)) ;    bp += (size_t)ne*nn*sizeof(double) ;

                eseed = 0 ;
                for(int i=0; i<3; i++) nseed[i][0] = nseed[i][1] = 0 ;

                T = new double [ne*nn*nn*nn] ;
                pre = new double [ne*nn*nn*nn] ;
                mub = new double [ne*nn*nn*nn] ;
                muq = new double [ne*nn*nn*nn] ;
                mus = new double [ne*nn*nn*nn] ;
                memcpy(pre, bp, npts*sizeof(double)) ; bp += npts*sizeof(double) ;
                memcpy(T,   bp, npts*sizeof(double)) ; bp += npts*sizeof(double) ;
                memcpy(mub, bp, npts*sizeof(double)) ; bp += npts*sizeof(double) ;
                memcpy(muq, bp, npts*sizeof(double)) ; bp += npts*sizeof(double) ;
                memcpy(mus, bp, npts*sizeof(double)) ;
                // same small-pressure clamp the ASCII parse loop applies
                for(size_t i=0; i<npts; i++)
                        if( pre[i] < 1e-18 ) pre[i] = 0. ;
                return ;
        }

        const char *sp = buf.data() ;
        const char *se = sp + buf.size() ;
        bool okread = true ;
        auto nextd = [&](double &v){
                while(sp<se && (*sp==' '||*sp=='\n'||*sp=='\t'||*sp=='\r'||*sp=='\v'||*sp=='\f')) ++sp ;
                std::from_chars_result r = std::from_chars(sp, se, v) ;
                if(r.ec != std::errc()){ v = 0 ; okread = false ; return ; }  // 0 on failure, like istream>>
                sp = r.ptr ;
        } ;
        auto nexti = [&](int &v){
                while(sp<se && (*sp==' '||*sp=='\n'||*sp=='\t'||*sp=='\r'||*sp=='\v'||*sp=='\f')) ++sp ;
                std::from_chars_result r = std::from_chars(sp, se, v) ;
                if(r.ec != std::errc()){ v = 0 ; okread = false ; return ; }
                sp = r.ptr ;
        } ;

        nextd(B) ; nextd(volex0) ; nextd(delta0) ; nextd(aaa) ; nextd(bbb) ;
        if( B > 0 )
        {
        if(fabs(B-_B)          >1e-5){ cout <<"\n\n\n      B =      "<<B     <<" instead of "<<_B     <<"; exiting\n\n\n\n" ; exit(0) ; }
        if(fabs(volex0-_volex0)>1e-5){ cout <<"\n\n\n      volex0 = "<<volex0<<" instead of "<<_volex0<<"; exiting\n\n\n\n" ; exit(0) ; }
        if(fabs(delta0-_delta0)>1e-5){ cout <<"\n\n\n      delta0 = "<<delta0<<" instead of "<<_delta0<<"; exiting\n\n\n\n" ; exit(0) ; }
        if(fabs(aaa-_aaa)      >1e-5){ cout <<"\n\n\n      aaa =    "<<aaa   <<" instead of "<<_aaa   <<"; exiting\n\n\n\n" ; exit(0) ; }
        if(fabs(bbb-_bbb)      >1e-5){ cout <<"\n\n\n      bbb =    "<<bbb   <<" instead of "<<_bbb   <<"; exiting\n\n\n\n" ; exit(0) ; }
        }
       
        nextd(emax) ; nextd(e0) ; nexti(ne) ; nexti(nn) ;

        egrid = new double [ne] ;
        ngrid = new double [ne*nn] ;

        for(int ixe=0;  ixe<ne;  ixe++)
                nextd(egrid[ixe]) ;

        for(int ixe=0;  ixe<ne;  ixe++)
        for(int ixnb=0; ixnb<nn; ixnb++)
                nextd(ngrid[ixe*nn+ixnb]) ;

        eseed = 0 ;
        for(int i=0; i<3; i++) nseed[i][0] = nseed[i][1] = 0 ;

        T = new double [ne*nn*nn*nn] ;
        pre = new double [ne*nn*nn*nn] ;
        mub = new double [ne*nn*nn*nn] ;
        muq = new double [ne*nn*nn*nn] ;
        mus = new double [ne*nn*nn*nn] ;

        int line=4;
        for(int ie=0; ie<ne; ie++)
        for(int inb=0; inb<nn; inb++)
        for(int inq=0; inq<nn; inq++)
        for(int ins=0; ins<nn; ins++){
                line++;
                const int idx = index(ie,inb,inq,ins) ;
                nextd(pre[idx]) ; nextd(T[idx]) ;
                nextd(mub[idx]) ; nextd(muq[idx]) ; nextd(mus[idx]) ;
                if(!okread){
                  cout<<"reading failed, line="<<line<<endl;
                  exit(13);
                }
                if( pre[idx] < 1e-18 ) pre[idx] = 0. ;
                //if(line==1000000)cout<<endl<<"EoS3f::EoS3f "<<line<<"  "<<index(ie,inb,inq,ins)<<"  "
                //                 << pre[index(ie,inb,inq,ins)] << "  " << T[index(ie,inb,inq,ins)] << "  " << mub[index(ie,inb,inq,ins)] << "  "
                //                 << muq[index(ie,inb,inq,ins)] << "  " << mus[index(ie,inb,inq,ins)] <<endl;
        }
}

EoS3f::~EoS3f()
{
        delete [] T ;
        delete [] pre ;
        delete [] mub ;
        delete [] muq ;
        delete [] mus ;
}

void EoS3f::eosranges(double &_emax, double &_e0, double &_nmax, double &_n0, int &_ne, int &_nn)
{
        _emax=emax;
        _e0=e0;
        _nmax=0; //not defined 
        _n0=0; //not defined
        _ne=ne;
        _nn=nn;
}

// Guarded walk from the previous interval instead of a full binary search.
// For a non-decreasing grid this yields exactly the interval the original
// binary search returned (max k in [0,ne-2] with egrid[k]<=e, clamped at the
// edges), for every input; successive calls query nearby e, so the walk is
// typically 0-2 steps.
void EoS3f::getue(double e, int &ixe, double &ue, int &iout){
        int k1, k2;
        iout=0;
        k1 = eseed ;
        for(;k1<ne-2 && e>=egrid[k1+1];) k1++ ;
        for(;k1>0 && e<egrid[k1];) k1-- ;
        eseed = k1 ;
        k2 = k1+1 ;
        ixe=k1;
        if(ixe>ne-2) /*{iout=1;return;}// */  ixe = ne - 2 ;
        double dxe = egrid[k2] - egrid[k1] ;
        double xem = e - egrid[k1] ;
        ue = xem/dxe ;
}

// Same guarded-walk replacement as getue(); rows may contain runs of equal
// values (e.g. all-zero rows at low e), for which the walk reproduces the
// binary search result (max k in [0,nn-2] with ngrid[row][k]<=n) exactly.
// islot selects the per-density seed (0=nb, 1=nq, 2=ns).
void EoS3f::getun(int ie, double n, int &ixn_low, int &ixn_up, double &un_low, double &un_up, int &iout, int islot){
        int k1, k2;
        const double *grlo = ngrid + ie*nn ;      // row ie
        const double *grup = grlo + nn ;          // row ie+1
        // --- lower side, so for ngird[ie][...]
        //if(n<grlo[0] )/*{iout=1;return;}*/ cout<<"getun " <<n<<"  "<<grlo[0]<<endl;
        //if(n>grlo[nn-1] )/*{iout=1;return;}*/ cout <<"getun " <<n<<"  "<<grlo[nn-1]<<endl;
        iout=0;
        if(n>grlo[nn-1]) n=grlo[nn-1];
        if(n<grlo[0]) n=grlo[0];
        k1 = nseed[islot][0] ;
        for(;k1<nn-2 && n>=grlo[k1+1];) k1++ ;
        for(;k1>0 && n<grlo[k1];) k1-- ;
        nseed[islot][0] = k1 ;
        k2 = k1+1 ;
        ixn_low=k1;
        if(ixn_low>nn-2) /*{iout=1;return;} // */ ixn_low = nn - 2 ;
        if(ixn_low<0) /*{iout=1;return;}  // */ ixn_low = 0 ;
        double dxn = grlo[k2] - grlo[k1] ;
        double xnm=n - grlo[k1] ;
        un_low = 0 ;
        if(dxn>0.)un_low = xnm/dxn ;
        // --- upper side so for ngird[ie+1][...]
        //if(n<grup[0] )/*{iout=1;return;}*/ cout <<"getun " <<n<<"  "<<grup[0]<<endl;
        //if(n>grup[nn-1] )/*{iout=1;return;}*/ cout <<"getun " <<n<<"  "<<grup[nn-1]<<endl;
        k1 = nseed[islot][1] ;
        for(;k1<nn-2 && n>=grup[k1+1];) k1++ ;
        for(;k1>0 && n<grup[k1];) k1-- ;
        nseed[islot][1] = k1 ;
        k2 = k1+1 ;
        ixn_up=k1;
        if(ixn_up>nn-2) /*{iout=1;return;}  // */ ixn_up = nn - 2 ;
        if(ixn_up<0) /*{iout=1;return;}  // */ ixn_up = 0 ;
        dxn = grup[k2] - grup[k1] ;
        xnm=n - grup[k1] ;
        un_up = 0 ;
        if(dxn>0.)un_up = xnm/dxn ;
        //cout<< "EoS3f::getun  "<<un_low<<" "<<un_up<<endl;
}

void EoS3f::eos(double e, double nb, double nq, double ns,
                double &_T, double &_mub, double &_muq, double &_mus, double &_p)
{
        if(e<0.) { _T = _mub = _muq = _mus = _p = 0. ; }
        int ixe, ixnb_low, ixnb_up, ixnq_low, ixnq_up, ixns_low, ixns_up, iout ;
        double ue, ub_low, ub_up, uq_low, uq_up, us_low, us_up ;
   
        getue(e,  ixe , ue, iout);                                 if(iout==1){_T = _mub = _muq = _mus = _p = 999.;return;}
        getun(ixe, nb, ixnb_low, ixnb_up, ub_low, ub_up, iout, 0); if(iout==1){_T = _mub = _muq = _mus = _p = 999.;return;}
        getun(ixe, nq, ixnq_low, ixnq_up, uq_low, uq_up, iout, 1); if(iout==1){_T = _mub = _muq = _mus = _p = 999.;return;}
        getun(ixe, ns, ixns_low, ixns_up, us_low, us_up, iout, 2); if(iout==1){_T = _mub = _muq = _mus = _p = 999.;return;}
                
        const double  we [2] = {1.-ue, ue} ;
        const double wnb [2][2] = {{1.-ub_low, ub_low},{1.-ub_up, ub_up}} ;
        const double wnq [2][2] = {{1.-uq_low, uq_low},{1.-uq_up, uq_up}} ;
        const double wns [2][2] = {{1.-us_low, us_low},{1.-us_up, us_up}} ;
        const int ixnb[2] = {ixnb_low, ixnb_up} ;
        const int ixnq[2] = {ixnq_low, ixnq_up} ;
        const int ixns[2] = {ixns_low, ixns_up} ;

        _T = _mub = _muq = _mus = _p = 0. ;
        for(int je=0; je<2; je++)
        for(int jnb=0; jnb<2; jnb++)
        for(int jnq=0; jnq<2; jnq++)
        for(int jns=0; jns<2; jns++){
          const int idx = index(ixe+je,ixnb[je]+jnb,ixnq[je]+jnq,ixns[je]+jns) ;
          if(pre[idx]==999.
             &&T[idx]==999.)goto NOSOL;
          const double w = we[je]*wnb[je][jnb]*wnq[je][jnq]*wns[je][jns] ;
          _p   += w*pre[idx] ;
          _T   += w*T[idx] ;
          _mub += w*mub[idx] ;
          _muq += w*muq[idx] ;
          _mus += w*mus[idx] ;
          // T 
          //cout << "EoS3f::eos  "<< ixe+je <<" "<<ixnb[je]+jnb <<" "<<ixnq[je]+jnq<<" "<<ixns[je]+jns  
          // <<"  " << index(ixe+je,ixnb[je]+jnb,ixnq[je]+jnq,ixns[je]+jns)<<"  " ;
          //cout << we[je]<<" "<<wnb[je][jnb]<<" "<<wnq[je][jnq]<<" "<<wns[je][jns]<<"        ";
          //cout << setw(10) << T[index(ixe+je,ixnb[je]+jnb,ixnq[je]+jnq,ixns[je]+jns)]<<"     ";
          //cout << egrid[ixe+je]<<" "<< ngrid[ixe+je][ixnb[je]+jnb] <<" "
          //<< ngrid[ixe+je][ixnq[je]+jnq]  <<" "<< ngrid[ixe+je][ixns[je]+jns]  <<endl;   
        }
        if(_p<0.) _p = 0. ;
        //cout << "EoS3f::eos "<<  e <<" "<< nb <<" "<< nq <<" "<< ns <<"     "<< _T<<endl;
        return;
        NOSOL: 
        //_T = _mub = _muq = _mus = 999.; _p = -999.;
        _T = _mub = _muq = _mus = _p = 0.;
}

// p-only variant of eos(): same lookups, same weights, same accumulation order
// for _p, so the returned value is bit-identical to eos()'s _p output; it only
// skips computing T/mub/muq/mus, which callers of p() discard.
double EoS3f::p(double e, double nb, double nq, double ns)
{
        int ixe, ixnb_low, ixnb_up, ixnq_low, ixnq_up, ixns_low, ixns_up, iout ;
        double ue, ub_low, ub_up, uq_low, uq_up, us_low, us_up ;

        getue(e,  ixe , ue, iout);                                 if(iout==1)return 999.;
        getun(ixe, nb, ixnb_low, ixnb_up, ub_low, ub_up, iout, 0); if(iout==1)return 999.;
        getun(ixe, nq, ixnq_low, ixnq_up, uq_low, uq_up, iout, 1); if(iout==1)return 999.;
        getun(ixe, ns, ixns_low, ixns_up, us_low, us_up, iout, 2); if(iout==1)return 999.;

        const double  we [2] = {1.-ue, ue} ;
        const double wnb [2][2] = {{1.-ub_low, ub_low},{1.-ub_up, ub_up}} ;
        const double wnq [2][2] = {{1.-uq_low, uq_low},{1.-uq_up, uq_up}} ;
        const double wns [2][2] = {{1.-us_low, us_low},{1.-us_up, us_up}} ;
        const int ixnb[2] = {ixnb_low, ixnb_up} ;
        const int ixnq[2] = {ixnq_low, ixnq_up} ;
        const int ixns[2] = {ixns_low, ixns_up} ;

        double _p = 0. ;
        for(int je=0; je<2; je++)
        for(int jnb=0; jnb<2; jnb++)
        for(int jnq=0; jnq<2; jnq++)
        for(int jns=0; jns<2; jns++){
          const int idx = index(ixe+je,ixnb[je]+jnb,ixnq[je]+jnq,ixns[je]+jns) ;
          if(pre[idx]==999.
             &&T[idx]==999.) return 0. ;
          const double w = we[je]*wnb[je][jnb]*wnq[je][jnq]*wns[je][jns] ;
          _p += w*pre[idx] ;
        }
        if(_p<0.) _p = 0. ;
        return _p ;
}

void EoS3f::eosorginal(double T, double mu_b, double mu_q, double mu_s, double &e, double& n_b, double& n_q, double& n_s, double& p) {
    e =0; n_b=0; n_q=0; n_s=0; p = 0.0;
    mix(T, mu_b, mu_q, mu_s, e, n_b, n_q, n_s, p) ;
    
}