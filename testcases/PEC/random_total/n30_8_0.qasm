OPENQASM 2.0;
include "qelib1.inc";
qreg q[30];
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
h q[20];
h q[21];
h q[22];
h q[23];
h q[24];
h q[25];
h q[26];
h q[27];
h q[28];
h q[29];
t q[13];
ccx q[9], q[20], q[6];
cx q[14], q[5];
ccx q[18], q[9], q[16];
ccx q[15], q[25], q[11];
h q[26];
h q[10];
t q[14];
s q[17];
t q[1];
h q[14];
s q[25];
cx q[20], q[28];
t q[20];
ccx q[4], q[5], q[18];
h q[21];
s q[29];
ccx q[21], q[22], q[5];
cx q[1], q[3];
ccx q[11], q[24], q[10];
cx q[11], q[7];
h q[26];
cx q[25], q[23];
cx q[23], q[19];
ccx q[16], q[28], q[2];
ccx q[19], q[22], q[18];
cx q[9], q[0];
t q[17];
cx q[3], q[21];
ccx q[23], q[9], q[27];
h q[19];
t q[1];
ccx q[24], q[5], q[25];
t q[21];
s q[10];
cx q[7], q[12];
t q[20];
ccx q[26], q[0], q[11];
t q[18];
cx q[6], q[14];
cx q[18], q[23];
ccx q[10], q[13], q[15];
s q[5];
cx q[1], q[24];
s q[1];
cx q[27], q[0];
s q[3];
s q[16];
s q[21];
h q[21];
h q[27];
cx q[6], q[16];
h q[5];
cx q[5], q[3];
cx q[27], q[10];
ccx q[24], q[7], q[21];
cx q[1], q[3];
h q[6];
ccx q[3], q[17], q[22];
h q[19];
h q[0];
s q[8];
cx q[18], q[16];
s q[9];
t q[26];
ccx q[18], q[19], q[7];
t q[8];
s q[29];
cx q[5], q[10];
s q[7];
s q[23];
ccx q[29], q[12], q[15];
cx q[1], q[18];
s q[19];
h q[9];
h q[26];
h q[13];
h q[4];
s q[29];
t q[2];
s q[0];
ccx q[14], q[12], q[13];
cx q[6], q[8];
h q[22];
t q[15];
cx q[10], q[19];
h q[21];
h q[27];
h q[18];
ccx q[15], q[14], q[6];
