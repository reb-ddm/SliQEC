OPENQASM 2.0;
include "qelib1.inc";
qreg q[16];
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
h q[10];
h q[11];
h q[12];
h q[13];
h q[14];
cx q[4], q[7];
h q[3];
cx q[3], q[6];
tdg q[6];
cx q[12], q[6];
t q[6];
cx q[3], q[6];
tdg q[6];
cx q[12], q[6];
t q[6];
cx q[12], q[3];
tdg q[3];
cx q[12], q[3];
t q[12];
t q[3];
h q[3];
t q[3];
cx q[1], q[2];
h q[13];
cx q[5], q[1];
h q[11];
h q[13];
s q[12];
cx q[4], q[6];
s q[8];
cx q[1], q[14];
h q[5];
cx q[5], q[1];
tdg q[1];
cx q[2], q[1];
t q[1];
cx q[5], q[1];
tdg q[1];
cx q[2], q[1];
t q[1];
cx q[2], q[5];
tdg q[5];
cx q[2], q[5];
t q[2];
t q[5];
h q[5];
s q[9];
t q[14];
t q[3];
h q[6];
s q[12];
t q[7];
s q[6];
t q[1];
h q[2];
cx q[2], q[4];
tdg q[4];
cx q[9], q[4];
t q[4];
cx q[2], q[4];
tdg q[4];
cx q[9], q[4];
t q[4];
cx q[9], q[2];
tdg q[2];
cx q[9], q[2];
t q[9];
t q[2];
h q[2];
s q[13];
s q[12];
cx q[0], q[6];
t q[9];
h q[6];
t q[8];
h q[4];
cx q[4], q[0];
tdg q[0];
cx q[12], q[0];
t q[0];
cx q[4], q[0];
tdg q[0];
cx q[12], q[0];
t q[0];
cx q[12], q[4];
tdg q[4];
cx q[12], q[4];
t q[12];
t q[4];
h q[4];
s q[0];
h q[12];
cx q[12], q[8];
tdg q[8];
cx q[0], q[8];
t q[8];
cx q[12], q[8];
tdg q[8];
cx q[0], q[8];
t q[8];
cx q[0], q[12];
tdg q[12];
cx q[0], q[12];
t q[0];
t q[12];
h q[12];
t q[4];
h q[0];
cx q[0], q[14];
tdg q[14];
cx q[7], q[14];
t q[14];
cx q[0], q[14];
tdg q[14];
cx q[7], q[14];
t q[14];
cx q[7], q[0];
tdg q[0];
cx q[7], q[0];
t q[7];
t q[0];
h q[0];
h q[8];
t q[12];
h q[0];
cx q[0], q[14];
tdg q[14];
cx q[9], q[14];
t q[14];
cx q[0], q[14];
tdg q[14];
cx q[9], q[14];
t q[14];
cx q[9], q[0];
tdg q[0];
cx q[9], q[0];
t q[9];
t q[0];
h q[0];
t q[7];
cx q[2], q[13];
s q[0];
cx q[14], q[13];
h q[14];
t q[10];
cx q[14], q[2];
h q[1];
s q[13];
t q[0];
cx q[1], q[0];
s q[0];
cx q[3], q[2];
x q[4];
cx q[5], q[6];
s q[5];
t q[6];
cx q[5], q[6];
t q[8];
t q[12];
h q[10];
t q[10];
t q[11];
t q[14];
h q[11];
t q[12];
cx q[15], q[13];
