  function  [sd,sl] = integratedsimilarity2(lncSim,disSim,id,il,dis_GS_Sim,lnc_GS_Sim,x)   
 

sd = x(1).*id+x(2).*disSim+x(3).*dis_GS_Sim;

sl = x(4).*il+x(5).*lncSim+x(6).*lnc_GS_Sim;


 end