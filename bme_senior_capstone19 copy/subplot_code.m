% File Name: subplut_code.m
% Authors: Maya Laughton and Gonzalo Miranda
% Affiliation: Tufts University, Biomedical Engineering Dpt.
% Class: BME 7 & 8
% Date: Fall 2018, Spring 2019
% Project: Senior Capstone
% Study Title: Learning Mechanisms of the Pavlovian pathway in the
%              amygdala: a computational modeling and recording study
% Description: In Silico Fear Conditioning Model. Models the neural state 
%              of a virtual rodent undergoing a fear conditioning 
%              experimental procedure
% Associated files: 
%   senior_capstone2019.m - Initializes parameter values. Calls
%                         network_math_v1, consolidate.m, subplot_code.m          
%   network_math_v2.m - The calculations to determine neural state of brain
%                       under experimentation.
%   consolidate_v1.m - For memory consolidation. To test a post-conditioned
%                   rat. Increases LA-CEm connectivity. Set conductivity of
%                   cxt cells. 
savefig('Fsmrd.fig')

% Load saved figures
subplot1=hgload('Events.fig');
subplot2=hgload('Fsmrd.fig');
subplot3=hgload('Vp (red) & Vq (blue)');
subplot4=hgload('Ap(red),Aq(blue),Af(green)');

% Prepare subplots
figure

h(4)=subplot(4,1,1);
ylim([0 1]);
xlim([0 Npoints]);
ylabel('Firing Rate')

h(2)=subplot(4,1,2);
ylim([0 1]);
xlim([0 Npoints]);
ylabel('Freezing')

h(3)=subplot(4,1,3);
ylim([0 80]);
xlim([0 Npoints]);
ylabel('mV')

h(1)=subplot(4,1,4);
ylim([0 0.9]);
xlim([0 Npoints]);
ylabel('Events')
xlabel('Seconds')

% Paste figures on the subplots
copyobj(allchild(get(subplot1,'CurrentAxes')),h(1));

copyobj(allchild(get(subplot2,'CurrentAxes')),h(2));

copyobj(allchild(get(subplot3,'CurrentAxes')),h(3));

copyobj(allchild(get(subplot4,'CurrentAxes')),h(4));