OPENQASM 2.0;
include "qelib1.inc";
qreg q[20];
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
h q[15];
h q[16];
h q[17];
h q[18];
h q[19];
t q[13];
t q[0];
s q[14];
s q[1];
t q[5];
s q[4];
h q[7];
s q[17];
h q[2];
cx q[15], q[16];
h q[18];
h q[18];
cx q[13], q[18];
s q[4];
s q[18];
cx q[17], q[6];
h q[19];
ccx q[14], q[19], q[12];
t q[19];
cx q[2], q[6];
ccx q[1], q[10], q[2];
h q[10];
t q[8];
s q[2];
t q[1];
ccx q[18], q[13], q[7];
cx q[5], q[11];
t q[17];
ccx q[18], q[19], q[5];
ccx q[5], q[7], q[0];
s q[7];
ccx q[6], q[14], q[4];
t q[12];
cx q[0], q[18];
h q[17];
h q[16];
t q[8];
h q[4];
cx q[9], q[14];
h q[3];
s q[2];
ccx q[6], q[2], q[1];
h q[3];
h q[7];
h q[14];
h q[14];
t q[3];
t q[5];
s q[11];
cx q[1], q[0];
t q[6];
ccx q[15], q[16], q[14];
t q[6];
h q[4];
h q[9];
t q[16];
cx q[14], q[18];
cx q[6], q[9];
cx q[16], q[11];
ccx q[16], q[3], q[7];
sdg q[0];
cx q[1], q[0];
tdg q[0];
cx q[1], q[0];
tdg q[0];
cx q[2], q[3];
tdg q[3];
cx q[2], q[3];
tdg q[3];
tdg q[2];
y q[4];
cx q[5], q[6];
tdg q[5];
tdg q[7];
tdg q[8];
cx q[7], q[8];
tdg q[7];
s q[19];
s q[12];
cx q[18], q[14];
h q[11];
ccx q[10], q[13], q[15];
cx q[13], q[11];
cx q[12], q[19];
ccx q[13], q[14], q[16];
cx q[16], q[13];
h q[10];
