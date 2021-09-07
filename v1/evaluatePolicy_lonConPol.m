function action1 = evaluatePolicy(observation1)
%#codegen

% Reinforcement Learning Toolbox
% Generated on: 19-May-2021 11:45:41

action1 = localEvaluate(observation1);
end
%% Local Functions
function action1 = localEvaluate(observation1)
persistent policy
if isempty(policy)
	policy = coder.loadDeepLearningNetwork('agentData_lonConPol.mat','policy');
end
observation1 = reshape(observation1,[1 1 7]);
action1 = predict(policy,observation1);
end