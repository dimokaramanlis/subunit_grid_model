function gforder = orderGratingFlicker(Nbase, oseed, Npresent)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


rnd_val = ran1(oseed, Npresent);
gforder = floor(Nbase*rnd_val) + 1;


end

