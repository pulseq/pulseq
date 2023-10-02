function init(varargin)
% mrMusic.init : create frequency and note variables
% creates global Matlab variables named after notes and octaves, for
% example c1, e2, etc; five octaves are created. The middle C is c1, the 
% higher octaves are nubered 2 and 3; left of the first octave is a small
% (suffix b) andd further left is bid (suffix bb). Sufixes 'is' and 'es'
% are added for 'sharp' and 'flat', respectively, eg. cis1, mes2, etc.
% The symbol 'o' is a delay.
% These variables are then used to create a melody in a humal-like
% notation. The melody is created bar-by-bar for the three voices 
% corresponding to the three gradient channels, with the duration specified 
% by dividing the desired note by a number, e.g. a1/16 is middle A with the
% duration of 1/16.
%
% Internally this syntax is achieved by defining the note variables as 
% complex numbers, where the complex part encodes the duration and the
% ratio of the real to the complex parts the pitch.
%
% It is possible to create non-standard tuning by setting the global Matlab 
% variable 'a_init' to a custom value (other than 440 Hz) prior to calling 
% this script for avoiding mechanicaal resonances for specific scanners
%

octave  ={'c','cis','d','dis','e','f','fis','g','gis','a','ais','h'};
synonyms={'cis','des'; 
          'dis','ees';
          'e',  'fes';
          'f',  'eis';
          'fis','ges';
          'ais','hes';}; % we cannot define sysnonyms for his and ces because this will break the octave limits
mult={'bb',0.25;'b',0.5;'1',1;'2',2;'3',4};

persistent parser
validTemperaments={'equal', 'meantone', 'meantone/4', 'Bach', 'BachWellTempered','Bach-Lehman'};

if isempty(parser)
    parser = inputParser;
    parser.FunctionName = 'mrMusic.init';
    parser.addParamValue('temperament',validTemperaments{1},...
        @(x) any(validatestring(x,validTemperaments)));
    parser.addParamValue('a_init',440,@isnumeric);
    %parser.addParamValue('setAsDefault',false,@islogical);
end
parse(parser,varargin{:});
opt = parser.Results;

htone=2^(1/12);
c_init=opt.a_init/htone^9;

temp_equal_cents=(0:11)*100;

% tuting tables. source: http://www.instrument-tuner.com/temperaments.html
delta_meantone=[8.798	-9.775	2.933	15.640	-2.932	11.731	-7.819	5.865	-10.752	0.000	13.686	-5.864];
delta_meantone4=[10.265	-13.686	3.422	-20.529	-3.421	13.686	-10.265	6.843	-17.108	0.000	-23.951	-6.843];
delta_bach_kellner_wt=[9.774	-0.001	3.258	3.909	-3.258	7.819	-1.956	6.516	1.954	0.000	5.864	-1.303];
delta_bach_kellner=[8.211	-1.564	2.737	2.346	-2.737	6.256	-3.519	5.474	0.391	0.000	4.301	-0.782];
delta_bach_lehman=[+5.9 +3.9 +2 +3.9 -2 +7.8 +2 +3.9 +3.9 0 +3.9 0];

switch opt.temperament 
    case 'equal'
        delta=zeros(1,12);
    case 'meantone'
        delta=delta_meantone;
    case 'meantone/4'
        delta=delta_meantone4;
    case 'Bach'
        delta=delta_bach_kellner;
    case 'BachWellTempered'
        delta=delta_bach_kellner_wt;
    case 'Bach-Lehman'
        delta=delta_bach_lehman;
end

for i=1:12
    %f=c_init*htone^(i-1);
    f=c_init*(2^((temp_equal_cents(i)+delta(i))/1200));
    for j=1:size(mult,1)
        fm=f*mult{j,2};
        eval([octave{i} mult{j,1} '=' num2str(fm) '+1j;']); % we add a complex unity to encode the note duration and allow the human-readable notation (see below)
        assignin('caller',[octave{i} mult{j,1}], fm+1j);
    end
end

% % Werckmeister (1691) / Bach (1722) wohltemperirt, source: http://plaza.ufl.edu/wnb/baroque_temperament.htm
% Werckmeister=[0.0 90.2 194.6 294.1 389.1 498.0 588.3 697.3 792.2 891.8 996.1 1091.1];
% delta_Werckmeister=Werckmeister(10)-Werckmeister-(9:-1:-2)*100;
% for i=1:12
%     f=c_init*(2^(Werckmeister(i)/1200));
%     for j=1:size(mult,1)
%         eval([octave{i} mult{j,1} '=' num2str(f*mult{j,2}) '+1j;']); % we add a complexunity to encode the note duration and allow the human-readable notation (see below)
%     end
% end

% % meantone table, source: Michtgan Tech University 
%meantone=[32.54	34.73	36.64	38.71	41.28	43.36	46.52	48.83	51.85	55	57.79	61.98	65.07	69.46	73.27	77.42	82.57	86.72	93.05	97.65	103.7	110	115.58	123.96	130.14	138.92	146.54	154.83	165.14	173.45	186.1	195.3	207.41	220	231.16	247.92	260.29	277.84	293.08	309.66	330.28	346.9	372.19	390.61	414.82	440	462.33	495.84	520.58	555.69	586.17	619.33	660.56	693.8	744.39	781.21	829.64	880	924.65	991.68	1041.16	1111.37	1172.34	1238.65	1321.12	1387.6	1488.78	1562.43	1659.28	1760	1849.31	1983.36	2082.31	2222.75	2344.68	2477.3	2642.24	2775.19	2977.56	3124.86	3318.56	3520	3698.61	3966.72	4164.63	4445.49	4689.36	4954.61];
% for j=1:size(mult,1)
%     for i=1:12
%         eval([octave{i} mult{j,1} '=' num2str(meantone(j*12+i)) '+1j;']); 
%     end
% end

for i=1:size(synonyms,1)
    for j=1:size(mult,1)
        %eval([synonyms{i,2}  mult{j,1} '=' synonyms{i,1}  mult{j,1} ';']);
        assignin('caller', [synonyms{i,2}  mult{j,1}], eval([synonyms{i,1}  mult{j,1}]))
    end
end

assignin('caller','o',0+1j); % delay (zero frequency)

% done with preparations! ready to define the melody! 
