function action1 = evaluatePolicy(observation1)
%#codegen

% Reinforcement Learning Toolbox
% Generated on: 22-Apr-2021 12:51:07

action1 = localEvaluate(observation1);
end
%% Local Functions
function action1 = localEvaluate(observation1)
persistent policy
if isempty(policy)
	policy = coder.loadDeepLearningNetwork('agentData_latConPol.mat','policy');
end
observation1 = reshape(observation1,[1 1 30]);
action1 = predict(policy,observation1);
end