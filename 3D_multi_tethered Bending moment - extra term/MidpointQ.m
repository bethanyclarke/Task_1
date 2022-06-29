function [out] = MidpointQ(q1,q2)
% Given 2 unit quaternions, q1 and q2, this function outputs the quaternion
% 'half way' between them in a rotational sense.

% First, we calculate the rotation mapping q1 to q2.

q1_conj = [q1(1),-q1(2),-q1(3),-q1(4)]; % Conjugate of q1.

q_rot = QuaternionProduct(q2,q1_conj); % q_rot is the quaternion mapping q1 to q2.

half_rot = unit_q_sqrt(q_rot); % Produce the half rotation.

out = QuaternionProduct(half_rot,q1); % Produce the midpoint.


end

