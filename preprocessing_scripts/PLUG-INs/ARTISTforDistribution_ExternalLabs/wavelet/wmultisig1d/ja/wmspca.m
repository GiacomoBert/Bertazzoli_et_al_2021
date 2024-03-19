% WMSPCA  �}���`�X�P�[���听������
%
%   [X_SIM,QUAL,NPC,DEC_SIM,PCA_Params] = WMSPCA(X,LEVEL,WNAME,NPC) �܂���
%   [...] = WMSPCA(X,LEVEL,WNAME,'mode',EXTMODE,NPC) �́A
%   �E�F�[�u���b�g�x�[�X�̃}���`�X�P�[�� PCA ���瓾������͍s�� X ��
%   �ȗ������ꂽ�s�� X_SIM ��Ԃ��܂��B
%
%   ���͍s�� X �́A������Ɋi�[���ꂽ���� N �̐M�� P ���܂�ł��܂��B
%   ������ N > P �ł��B
%
%   �E�F�[�u���b�g�����p�����[�^
%   ----------------------------
%   �E�F�[�u���b�g�����́A�������x�� LEVEL �ƃE�F�[�u���b�g WNAME ���g����
%   �s�Ȃ��܂��BEXTMODE �́ADWT �ɑ΂���g�����[�h�ł� (�f�t�H���g�� 
%   DWTMODE �ŕԂ���郂�[�h�ł�)�B
%
%   MDWTDEC ���g���ē����镪�� DEC �����p�\�ȏꍇ�A
%   [...] = WMSPCA(X,LEVEL,WNAME,'mode',EXTMODE,NPC)
%   �̑���� [...] = WMSPCA(DEC,NPC) ���g�p���邱�Ƃ��\�ł��B
%
%   �听���p�����[�^ NPC
%   --------------------
%   NPC ���x�N�g���̏ꍇ�ALEVEL+2 �̒����łȂ���΂Ȃ�܂���B���s�����
%   �e PCA �ɑ΂���c��̎听���̐����܂݂܂��B
%      - NPC(d) �́A1 <= d <= LEVEL �ɑ΂��ă��x�� d �ɂ����� detail ��
%        �c��̔񒆐S���̎听�����ł��B
%      - NPC(LEVEL+1) �́A���x�� LEVEL �ɂ����� approximations �ɑ΂���
%        �c��̔񒆐S���̎听�����ł��B
%      - NPC(LEVEL+2) �́A�E�F�[�u���b�g�č\����� PCA �ɑ΂���c���
%        �听�����ł��B
%      NPC �́A1 <= d <= LEVEL+2 �ɑ΂��� 0 <= NPC(d) <= P �łȂ����
%      �Ȃ�܂���B
%
%   NPC = 'kais' (�܂��� 'heur') �̏ꍇ�A�c��̎听�����́AKaiser ��
%   ���[�� (�܂��͌o����) ���g���Ď����I�ɑI������܂��B
%      - Kaiser �̃��[���́A���ׂĂ̌ŗL�l�̕��ς𒴂���ŗL�l�ɑ΂���
%        ������ێ����܂��B
%      - �o�����́A���ׂĂ̌ŗL�l�̘a�� 0.05 �{�𒴂���ŗL�l�ɑ΂���
%        ������ێ����܂��B
%   NPC = 'nodet' �̏ꍇ�Adetail �� "�폜" ����Aapproximation ���c��܂��B
%
%   �o��
%   ----
%   X_SIM �́A�s�� X �̊ȗ������ꂽ�s��ł��B
%   QUAL �́A�p�[�Z���g�ŕ\�L�ő��Ε��ϓ��덷�ŗ^������č\���M����
%   �i�����܂ޒ��� P �̃x�N�g���ł��B
%   NPC �́A�I�������c��̎听�����̃x�N�g���ł��B
%   DEC_SIM �́AX_SIM �̃E�F�[�u���b�g�����ł��B�����̍\���̂̂��ڂ���
%   ���ɂ��Ă� MDWTDEC ���Q�Ƃ��Ă��������B
%   PCA_Params �́A���̂悤�� LEVEL+2 �̒����̍\���̔z��ł��B
%     PCA_Params(d).pc = PC�B�����ŁA
%        PC �́A�听���� P�~P �̍s��ł��B
%        ��́A���U�̍~���ɏ]���Ċi�[���ꂽ���̂ł��B
%     PCA_Params(d).variances = VAR�B�����ŁA
%        VAR �́A�听�����͂̕��U�x�N�g���ł��B
%     PCA_Params(d).npc = NPC.


%   Copyright 1995-2007 The MathWorks, Inc.
