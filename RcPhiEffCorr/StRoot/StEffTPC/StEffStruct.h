struct McEvent
{
  float runId;
  float eventId;
  float mcVx;
  float mcVy;
  float mcVz;
  float vx;
  float vy;
  float vz;
  float vzVpd;
  float centrality;
  float gRefMult;
  float RefMult;
  float posRefMult;
  float negRefMult;
  float zdc;
  float bbc;
  float nMcTracks;
  float nRTracks;
  float magField;
  float t0;
  float t1;
  float t2;
  float t3;
  float t4;
  float t5;
};

struct McVecMeson
{
  float pt;
  float p;
  float eta;
  float y;
  float phi;
  float label;
  float geantId;
};

struct McDecayDau
{
  float pt;
  float eta;
  float phi;
  float geantId;
  float startVtxX;
  float startVtxY;
  float startVtxZ;
  float stopVtxX;
  float stopVtxY;
  float stopVtxZ;
};

struct RcDecayDau
{
  float pt;
  float eta;
  float phi;
  float nFit;
  float nMax;
  float nCom;
  float nDedx;
  float Dedx;
  float nSigKP;
  float nSigKM;
  float dca;
  float dcaXY;
  float dcaZ;
  float McPt;
  float McEta;
  float McPhi;
};

struct RcVecMeson
{
  float InvMass;
  float pt;
  float eta;
  float y;
  float phi;
  float McPt;
  float McEta;
  float McPhi;
};
