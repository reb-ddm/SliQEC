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
s q[11];
s q[9];
s q[1];
ccx q[13], q[28], q[17];
t q[6];
s q[9];
h q[4];
h q[10];
h q[4];
h q[16];
h q[3];
h q[3];
s q[12];
ccx q[9], q[0], q[14];
h q[6];
ccx q[23], q[19], q[16];
ccx q[0], q[25], q[7];
s q[0];
ccx q[14], q[8], q[27];
t q[28];
ccx q[19], q[26], q[5];
ccx q[19], q[12], q[28];
t q[3];
s q[3];
ccx q[4], q[15], q[0];
ccx q[1], q[9], q[6];
h q[8];
t q[14];
cx q[28], q[25];
s q[8];
ccx q[17], q[13], q[28];
cx q[9], q[14];
s q[19];
ccx q[13], q[3], q[12];
s q[3];
t q[26];
ccx q[7], q[28], q[21];
h q[16];
h q[0];
cx q[11], q[6];
s q[10];
cx q[3], q[1];
cx q[13], q[27];
t q[8];
cx q[19], q[8];
cx q[18], q[4];
s q[0];
h q[5];
h q[10];
ccx q[11], q[16], q[14];
cx q[25], q[6];
cx q[1], q[11];
ccx q[12], q[21], q[14];
ccx q[19], q[12], q[18];
t q[26];
ccx q[18], q[3], q[20];
s q[9];
s q[20];
ccx q[24], q[18], q[23];
ccx q[12], q[8], q[4];
ccx q[2], q[28], q[9];
cx q[9], q[24];
cx q[19], q[0];
cx q[6], q[16];
t q[25];
s q[0];
ccx q[21], q[16], q[9];
h q[13];
t q[23];
cx q[29], q[11];
ccx q[23], q[26], q[12];
s q[10];
h q[1];
h q[19];
ccx q[22], q[20], q[13];
h q[13];
t q[20];
s q[10];
t q[22];
t q[11];
s q[10];
t q[2];
t q[10];
t q[20];
t q[28];
s q[13];
ccx q[12], q[7], q[19];
h q[7];
t q[14];
s q[25];
