struct McEvent
{
  float runId;
  float eventId;
  float McVx;
  float McVy;
  float McVz;
  float vx;
  float vy;
  float vz;
  float vzVpd;
  float Centrality;
  float gRefMult;
  float RefMult;
  float posRefMult;
  float negRefMult;
  float zdc;
  float bbc;
  float nMcTracks;
  float nRcTracks;
  float magField;
  float t0;
  float t1;
  float t2;
  float t3;
  float t4;
  float t5;
};

struct McDecayDau
{
  float McPt;
  float McP;
  float McEta;
  float McY;
  float McPhi;
  float McGeantId;
  float eventGenLabel;
  float StartVtxX;
  float StartVtxY;
  float StartVtxZ;
  float StopVtxX;
  float StopVtxY;
  float StopVtxZ;
};

struct RcDecayDau
{
  float gRcPt; // global RcTracks
  float gRcEta;
  float gRcPhi;
  float gRcNfit;
  float gRcNmax;
  float gRcNcom;
  float gRcNdedx;
  float gRcDedx;
  float gRcNsigKP;
  float gRcNsigKM;
  float gRcDca;
  float gRcDcaXY;
  float gRcDcaZ;

  float pRcPt; // primary RcTracks
  float pRcEta;
  float pRcPhi;
  float pRcNfit;
  float pRcNmax;
  float pRcNcom;
  float pRcNdedx;
  float pRcDedx;
  float pRcNsigKP;
  float pRcNsigKM;
  float pRcDca;
  float pRcDcaXY;
  float pRcDcaZ;

  float McPt;
  float McEta;
  float McPhi;
};
