OPENQASM 2.0;
include "qelib1.inc";
qreg q[11];
h q[0];
h q[1];
h q[2];
h q[3];
h q[4];
h q[5];
h q[6];
h q[7];
h q[8];
h q[9];
h q[1];
cx q[1], q[4];
tdg q[4];
cx q[7], q[4];
t q[4];
cx q[1], q[4];
tdg q[4];
cx q[7], q[4];
t q[4];
cx q[7], q[1];
tdg q[1];
cx q[7], q[1];
t q[7];
t q[1];
h q[1];
h q[7];
h q[4];
cx q[4], q[5];
tdg q[5];
cx q[6], q[5];
t q[5];
cx q[4], q[5];
tdg q[5];
cx q[6], q[5];
t q[5];
cx q[6], q[4];
tdg q[4];
cx q[6], q[4];
t q[6];
t q[4];
h q[4];
h q[6];
h q[8];
cx q[5], q[0];
h q[1];
cx q[1], q[5];
tdg q[5];
cx q[7], q[5];
t q[5];
cx q[1], q[5];
tdg q[5];
cx q[7], q[5];
t q[5];
cx q[7], q[1];
tdg q[1];
cx q[7], q[1];
t q[7];
t q[1];
h q[1];
h q[4];
cx q[4], q[0];
tdg q[0];
cx q[8], q[0];
t q[0];
cx q[4], q[0];
tdg q[0];
cx q[8], q[0];
t q[0];
cx q[8], q[4];
tdg q[4];
cx q[8], q[4];
t q[8];
t q[4];
h q[4];
h q[7];
cx q[7], q[8];
tdg q[8];
cx q[1], q[8];
t q[8];
cx q[7], q[8];
tdg q[8];
cx q[1], q[8];
t q[8];
cx q[1], q[7];
tdg q[7];
cx q[1], q[7];
t q[1];
t q[7];
h q[7];
h q[8];
h q[5];
cx q[5], q[9];
tdg q[9];
cx q[3], q[9];
t q[9];
cx q[5], q[9];
tdg q[9];
cx q[3], q[9];
t q[9];
cx q[3], q[5];
tdg q[5];
cx q[3], q[5];
t q[3];
t q[5];
h q[5];
cx q[8], q[0];
h q[0];
t q[6];
cx q[6], q[8];
s q[4];
h q[3];
s q[7];
t q[7];
s q[6];
t q[1];
s q[0];
h q[3];
cx q[3], q[7];
tdg q[7];
cx q[0], q[7];
t q[7];
cx q[3], q[7];
tdg q[7];
cx q[0], q[7];
t q[7];
cx q[0], q[3];
tdg q[3];
cx q[0], q[3];
t q[0];
t q[3];
h q[3];
h q[5];
s q[2];
cx q[4], q[8];
t q[0];
h q[3];
cx q[6], q[7];
h q[6];
cx q[6], q[0];
tdg q[0];
cx q[7], q[0];
t q[0];
cx q[6], q[0];
tdg q[0];
cx q[7], q[0];
t q[0];
cx q[7], q[6];
tdg q[6];
cx q[7], q[6];
t q[7];
t q[6];
h q[6];
t q[4];
t q[3];
cx q[3], q[4];
h q[9];
t q[8];
h q[5];
cx q[5], q[7];
tdg q[7];
cx q[6], q[7];
t q[7];
cx q[5], q[7];
tdg q[7];
cx q[6], q[7];
t q[7];
cx q[6], q[5];
tdg q[5];
cx q[6], q[5];
t q[6];
t q[5];
h q[5];
t q[6];
s q[9];
cx q[10], q[0];
