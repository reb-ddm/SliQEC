OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
t q[2];
cx q[2], q[4];
cx q[2], q[0];
cx q[2], q[0];
cx q[1], q[2];
s q[1];
t q[3];
s q[0];
cx q[3], q[4];
ccx q[3], q[1], q[2];
cx q[4], q[2];
cx q[1], q[2];
ccx q[0], q[4], q[1];
t q[3];
s q[3];
sdg q[0];
cx q[1], q[0];
sdg q[1];
sdg q[0];
h q[3];
s q[4];
t q[2];
cx q[5], q[4];
