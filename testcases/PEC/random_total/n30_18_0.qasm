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
cx q[11], q[15];
ccx q[13], q[8], q[5];
h q[29];
cx q[8], q[1];
h q[19];
h q[25];
t q[26];
ccx q[29], q[15], q[27];
cx q[26], q[29];
ccx q[28], q[20], q[15];
h q[6];
h q[11];
cx q[28], q[8];
t q[16];
ccx q[10], q[26], q[4];
t q[4];
h q[23];
t q[2];
t q[29];
t q[12];
h q[27];
t q[18];
s q[26];
cx q[16], q[9];
t q[11];
s q[18];
h q[3];
s q[15];
h q[10];
s q[8];
cx q[5], q[13];
s q[17];
ccx q[17], q[14], q[5];
cx q[22], q[4];
cx q[15], q[17];
t q[21];
ccx q[21], q[16], q[24];
h q[6];
s q[2];
t q[12];
cx q[12], q[19];
ccx q[12], q[1], q[18];
t q[4];
s q[8];
s q[18];
s q[27];
cx q[24], q[26];
h q[19];
t q[26];
ccx q[26], q[28], q[19];
s q[17];
ccx q[12], q[10], q[4];
ccx q[16], q[25], q[19];
cx q[16], q[3];
cx q[26], q[12];
t q[10];
cx q[18], q[21];
t q[12];
cx q[3], q[27];
s q[15];
cx q[14], q[5];
s q[23];
s q[22];
h q[20];
t q[4];
ccx q[16], q[5], q[18];
ccx q[14], q[7], q[17];
cx q[15], q[20];
s q[1];
ccx q[16], q[8], q[9];
ccx q[3], q[25], q[0];
cx q[8], q[13];
h q[28];
t q[15];
h q[26];
ccx q[26], q[15], q[24];
cx q[18], q[16];
t q[20];
cx q[13], q[17];
t q[19];
h q[28];
h q[2];
s q[6];
s q[8];
t q[20];
ccx q[27], q[5], q[17];
cx q[7], q[16];
h q[0];
t q[0];
s q[9];
