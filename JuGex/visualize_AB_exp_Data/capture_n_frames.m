function [F,frame_counter] = capture_n_frames( n,frame_counter,F,h2capture)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    for i=1:n
        F(frame_counter) = getframe(h2capture);
        frame_counter=frame_counter+1;
    end
end

