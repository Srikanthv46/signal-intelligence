import numpy as np
from spectrum_intel.receivers.rc_link import synth as S
from ber import OFS, demod_disc_int, SPS, _kmeans2
from coherent import detect_coh3 as detect_coh2, TAB
print(np.degrees(TAB).round(1))
def chan(snr,rng,dev,drift_hz_per_us,nbits=400):
    aa=[(0x8E89BED6>>i)&1 for i in range(32)]; pre=[aa[0]^(i%2==0) for i in range(8)]
    bits=np.r_[pre,aa,rng.integers(0,2,nbits-40)].astype(np.uint8)
    iq,_=S.gfsk_burst(1e6,dev,bits=bits,bt=0.5,cfo_hz=rng.uniform(-40e3,40e3),fs=OFS,rng=rng)
    t=np.arange(len(iq))/OFS
    iq=iq*np.exp(1j*np.pi*drift_hz_per_us*1e6*t**2)
    n0=1/(1e6*10**(snr/10))
    w=(rng.standard_normal(len(iq))+1j*rng.standard_normal(len(iq)))*np.sqrt(n0*OFS/2)
    Y=np.fft.fft(iq+w); f=np.fft.fftfreq(len(iq),1/OFS); Y[np.abs(f)>0.75e6]=0
    return np.fft.ifft(Y),bits
for dev,dr in ((250e3,0),(225e3,0),(275e3,0),(250e3,100),(250e3,400),(225e3,400)):
  for snr in (6,8,10):
    rng=np.random.default_rng(1); e=ed=n=0; hs=[]
    for t in range(40):
        y,bits=chan(snr,rng,dev,dr)
        a,_,h,_=detect_coh2(y,-9,len(bits),16); hs.append(h)
        tb=bits[1:]; m=len(tb)-4; e+=np.sum(a[40:m]!=tb[40:m]); n+=m-40
        v=demod_disc_int(y); idx=np.arange(len(bits))*SPS+4; idx=idx[idx<len(v)]; vv=v[idx]
        ed+=np.sum((vv>np.mean(_kmeans2(vv)))[41:m+1]!=bits[41:m+1].astype(bool))
    print(f"dev {dev/1e3:.0f}k drift {dr} Hz/us snr {snr}: coherent {e/n:.5f}  diff {ed/n:.5f}  h {np.mean(hs):.3f}")
