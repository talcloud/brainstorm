%function tree=createHeaderLena(xml_file,data_format,data_type,data_size, data_blocks, sensorsList, coilList, correctionList,sensorSamples,frequenciesList, frequencySamples,time_samples,sample_rate, pre_trigger )
%
 

function createFolderLENAHead(filename,header, overwrite )




 
% % % 
% %
% %
% % This function will create a LENA xml file with minimal mandatory
% % arguments
% %
% %   Inputs :
% %       - xml_file      ( string )  : file name to create
% %
% %       - data_format   ( string )  : endianness ('LittleEndian' or 'BigEndian')
% %
% %       - data_type     ( string )  : data type ( 'float', 'int', etc ... )
% 
% %       - data_size     ( string )  : size of data
% 
% %       - data_blocks               :  Cell array of string defining class name for each block/trial in the buffer 
% %                                       example data_blocks={'trial 1'; 'trial 2'; 'trial 3'; 'trial 4'};
% 
% %       - sensorList                : Struct made of 2 fields
% %                                           1. name:        Cell array of strings defining the names of the sensors
% %                                           2. category:    Cell array of strings specifying type of sensors (meg, eeg, adc, eog, ecg)
% 
% %                                   example:    sensorsList.name= {'PO7'; 'PO6';'PO5'};
% %                                               sensorsList.category= {'eeg'; 'meg' ;'eog'};
% 
% %       - coilList                  : struct  made of 3 fields: 
% 
% %                                       1. index:array of integer, the ith element of the array is the index of the sensor to which belongs the i th coil   (for MEG
% %                                       several coils can be specified)
% 
% %                                       2. position: array of float vertors, the ith vector is the position of the ith coil 
% 
% %                                       3. orientation: array of float vertors, the ith vector is the orientation of the ith coil 
% %
% %                                       example:    coilList.index=[2,2,3];
% %                                                   coilList.position={[0.23,1.85,0.36], [1.54,0.96,1.42], [1.87,1.69,1.32]};
% %                                                   coilList.orientation={[0,1,1], [1,0,1], [0,1,0]};
% %
% %       - correctionList            : struct made of 3 fields
% %               
% %                                      1. type: cell array of strings  specifying the type of correction
% 
% %                                      2. index: array of integer, the ith
% %                                      element of the array is the index of
% %                                      the sensor to which belongs the i th correction 
% %                                       (several corrections can be specified)
% %
% %                                      3. Element: Cell aray of element.
% %                                      Each element can be a singleton or
% %                                      an array of multiple singletons. 
% %                                      A singleton is a cell of 2 values:
% %                                      a string (name of the sensor) and a float (
% %                                      correction coefficient ) example: {'BG1-610',0.1328}
% %                   
% %                                       example:correctionList.type={'RB1G', 'IO3G', 'IO2G'};
% 
% %                                               correctionList.index=[1, 1, 3];
% 
% %                                               correctionList.Element={{'BG1-610',0.1328}; [{'BG2-610', 0.236}; {'BG3-610', 1.25}; {BG4-610', 4.53}]; {'BG5-610', 1.45}};
% % 
% % 
% % 
% 
% %       - sensorSamples             :   struct made of 4 fields
% %                                       1. list_sensors: cell array of cell
% %                                       elements. Each element is a
% %                                       supersensor
% %                                       
% %                                        2. unit: array of cell element: the
% %                                       first component is the unit of the
% %                                       sensor and the second is the index
% %                                       of the supersensor in the
% %                                       list_sensors
% %
% %                                        3.scale:array of cell element: the
% %                                       first component is the unit of the
% %                                       sensor and the second is the index
% %                                       of the super in the
% %                                       list_sensors
% %
% %                                       4. num_sensors:number of sensors in
% %                                       each supersensor of list_sensors
% %                                       
% %                                       example:
% %                                       sensorSamples.unit=[{'v',1};{'v',2}];
% %                                       sensorSamples.scale=[{0,1};{0,2}];
% %                                       sensorSamples.num_sensors=[1;2;1];
% %                                       sensorSamples.list_sensors={{'MZO02-610'};{'MZO01-610', 'MZO03-610'};{'MZO04-610'}};
% %
% %
% %       - frequenciesList           : Array of float defining frequencies
% %                                      example: frequenciesList=[10 20 30 40 50 60 70];
% %
% % 
% %       - frequencySamples          : Cell array of superfrequencies
% %                                       example: frequencySamples={{10, 20} ;{20} ;{30} ;{40, 50}; {50}; {60} ;{70}};
% %
% %       - timeSamples               : number of time samples (must be
% %                                         positive)
% %
% %       - sample_rate               : Sample rate of the time dimension
% %       - pre_trigger               : Pre trigger time (must be positive)
% %  
% %
% % Outputs :
% %       - tree          ( xmltree ) : xml object that has just been stored in the xml
% %       file.
% %
% % Notice : The parameters data_format, data_type, data_size are part of
% % the LENA data format. Please refer to its documentation for more
% % informations.
% %
% % See also : xmltree, insertSimpleElements, insertTimeRange,
% % copyLENAelement, deleteLENAelement
% %
% %
% 
% 
% %
% % high_filters: Array of float defining high pass filter values used on the
% %               the corresponding sensor
% %
% % low_filters: Array of float defining low pass filter values used on the
% %               the corresponding sensor
% % 
% % impedance: Array of float defining impedance of each sensor
% %
% % phys_max: Array of float defining maximum physical value for each sensor
% %
% % phys_min: Array of float defining minimum physical value for each sensor
% %
% 
% 
% 
% 
% 
% 
% % Author :
% %
% % Lydia Yahia
% % 05 03 2008
% %
% % LENA CNRS UPR640
% % 


%error(nargchk(3, 4, nargin))
n=nargin;
if (n==2)
    overwrite=false;
%else overwrite=true;
end
%create folder


if (~ isdir(filename)) || (isdir(filename) & overwrite)
    status=mkdir(filename);
    if (status )
    xml_file=fullfile(filename,'data.header');
    createHeader(xml_file,header );
    %WriteDataLENA(MatData, xml_file);
    else error('bad directory');
    return;
    end
else error('Already existing Directory')


%end
end
% error(nargchk(2,Inf,nargin,'struct'));
% 
% datablocks=     [];
% sensorList =    [];
% coilList=       []; 
% correctionList= []; 
% sensorSamples=  []; 
% frequenciesList=[];
% frequencySamples = [];
% time_samples=     [];
% sample_rate=        []; 
% pre_trigger=        [];
% 
% if ~isempty(header) 
%    [data_format,data_type,data_size, datablocks,sensorList,  sensorSamples, frequenciesList, frequencySamples ,time_samples,sample_rate, pre_trigger, history]= checkargt(header);
%     
% 
% end
% 
% 
% 
% % Create an xml tree
% tree = xmltree;
% 
% % Set the root element to 'LENA' :
% tree = set(tree,root(tree),'name','LENA');
% 
% % Add mandatory arguments :
% [ tree , format_uid ] = add(tree,root(tree),'element','data_format');
% tree = add(tree,format_uid,'chardata',data_format);
% 
% [ tree , type_uid ] = add(tree,root(tree),'element','data_type');
% 
% tree = add(tree,type_uid,'chardata',data_type);
% 
% [ tree , size_uid ] = add(tree,root(tree),'element','data_size');
% tree = add(tree,size_uid,'chardata',num2str(data_size));
% 
% if exist('history')
%     if ~(isempty(history))
%         tree=setHistory(tree, history);
%     end
% end
% 
% 
% if (exist('datablocks'))
% if ~(isempty(datablocks))
% tree=insertDataBlockRange(tree, datablocks);
% end;
% end;
% 
% 
% %add sensors and sensor samples
% if exist('sensorList') || exist('sensorsSamples')
% if (~isempty(sensorList) || ~isempty(sensorSamples))
%     if (  isfield(sensorList,'coils') &  ~isempty(sensorList.coils)) 
%         coilList=sensorList.coils;
%     end
% %     if (  isfield(sensorList,'category') &  ~isempty(sensorList.category)) 
% %         List=sensorList.coils;
% %     end
%     
%     
% 
% tree=insertSensorRange(tree,sensorList,sensorSamples, coilList, correctionList);
% end
% end
% 
% % add frequencies and frequencies samples
% % 
% 
% 
% if exist('frequenciesList') || exist('frequencySamples')
% if (~isempty(frequenciesList) || ~isempty(frequencySamples))
%     %frequencySamples=frequenciesList;
% tree = insertFrequencyRange(tree, frequenciesList, frequencySamples);
% 
% end
% end
% % Add timeRange
% if exist('time_samples') || exist('sample_rate') || exist('pre_trigger')
% if (~isempty(time_samples) & ~isempty(sample_rate) & ~isempty(pre_trigger))
% tree=insertTimeRange(tree,num2str(time_samples),num2str (sample_rate),num2str(pre_trigger));
% end
% end
% 
% 
% 
% 
% % Save xml :
% savexmlLena(tree,xml_file);
% %xmlwrite(xml_file,tree);
% 
% %return
% 
%   end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function  [data_format,data_type,data_size, datablocks,sensorList, sensorSamples, frequenciesList, frequencySamples ,time_samples,sample_rate, pre_trigger, history]= checkargt(header)
%     
% data_format=[];
% data_type=[];
% data_size=[]; 
% history=[];
% datablocks=[];    
% sensorList =[];
% coilList=[]; 
% correctionList=[]; 
% sensorSamples=[]; 
% frequenciesList=[];
% frequencySamples =[];
% time_samples= [];
% sample_rate= []; 
% pre_trigger=[];
%     
%     
%     i=1;
% %len=numel(argval{1});
% %while i<=  len
% %va1 = varargin{1,1}{i,1};
%  %  if ischar (va1) && ~isempty(va1)
%   %  a= strmatch(va1,{'sensorList', 'coilList', 'correctionList', 'sensorSamples', 'frequenciesList',' frequencySamples ','time_samples','sample_rate', 'pre_trigger' }, 'exact'))
%    % if ~isempty(a)
%     %    val=a[1];
%     
%     
%     if     isfield(header, 'datablocks')  & ~isempty(header.datablocks)
%               datablocks=header.datablocks;
%     end
% 
%     if   isfield(header, 'sensors')  &  ~ isempty(header.sensors)
%               sensorList=header.sensors;
%     end      
%     
%     if  isfield(header, 'sensorSamples')  &  ~ isempty(header.sensorSamples)
%               sensorSamples=header.sensorSamples;
%     end 
%     
%     
%      if  isfield(header, 'timeSamples')  &  ~ isempty(header.timeSamples)
%               time_samples=header.timeSamples;
%      end 
%     
%      
%      
%           if  isfield(header, 'sampleRate')  &  ~ isempty(header.sampleRate)
%               sample_rate=header.sampleRate;
%         end 
% 
%      if  isfield(header, 'preTrigger')  &  ~ isempty(header.preTrigger)
%               pre_trigger=header.preTrigger;
%      end 
%         
%      
%           if  isfield(header, 'frequencies')  &  ~ isempty(header.frequencies)
%               frequenciesList=header.frequencies;
%           end 
% 
%         
%           
%        if  isfield(header, 'frequencysamples')  &  ~ isempty(header.frequencysamples)
%               frequencySamples=header.frequencysamples;
%        end 
% 
% 
%       if  isfield(header, 'data_format')  &  ~ isempty(header.data_format)
%               data_format=header.data_format;
%       end 
%  
%       
%         
%       if  isfield(header, 'data_type')  &  ~ isempty(header.data_type)
%               data_type=header.data_type;
%       end 
%         
%       
%       if  isfield(header, 'data_size')  &  ~ isempty(header.data_size)
%               data_size=header.data_size;
%       end 
%  
%  
%         
%        if  isfield(header, 'history')  &  ~ isempty(header.history)
%               history=header.history;
%         end 
%  
%        
% %           case 'sensorList'
% %               sensorList=argval{1,1}{i+1,1};
% %               
% %            case 'coilList'
% %             coilList=argval{1,1}{i+1,1};
% %             
% %             case 'correctionList'
% %               correctionList=argval{1,1}{i+1,1};
% %               
% %                     case 'sensorSamples'
% %               sensorSamples=argval{1,1}{i+1,1};
% %               
% %                     case 'correctionList'
% %               correctionList=argval{1,1}{i+1,1};              
% %                     case 'frequenciesList'
% %               frequenciesList=argval{1,1}{i+1,1};
% %               
% %                     case 'frequencySamples'
% %               frequencySamples=argval{1,1}{i+1,1};
% %               
% %                     case 'time_sample'
% %               time_samples=argval{1,1}{i+1,1};
% %               
% %                     case 'sample_rate'
% %               sample_rate=argval{1,1}{i+1,1};
% %               
% %                     case 'pre_trigger'
% %               pre_trigger=argval{1,1}{i+1,1};
% %               
% %                  case 'history'
% %               history=argval{1,1}{i+1,1};
% %       end
%     %end
%  %  end
% 
% 
% %i=i+2;
% %end
%  %end
% end
