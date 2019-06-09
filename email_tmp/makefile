CXX      = g++
#CXXFLAGS = -std=c++11 -g -lm # 
#CXXFLAGS = -std=c++11 -g -lmat -l mx
#CXXFLAGS = -std=c++11 -g -lmat -lm
CXXFLAGS = -std=c++11 -g -leng -lmat -lmex -lut -lmx
INC_PATH = -I/home/abcd/Documents/LDPC_Density_Evolution/my_discreted_density_evolution/inc \
	       -I/usr/tools/matlab/R2017b/extern/include \
	   	   -I/usr/lib/Eigen/Eigen
LIBS     = -L/home/abcd/Documents/LDPC_Density_Evolution/my_discreted_density_evolution \
	       -L/usr/tools/matlab/R2017b/bin/glnxa64
SRCFILES = ./src/*.cpp                                		    # 指定目录 ./ 下的所有后缀是cpp的文件全部展开。
OBJS     = $(SRCFILES:.cpp=.o)                                  # OBJS将$(SRCS)下的.cpp文件转化为.o文件
TARGET   = integral                                             # 输出程序名称

all: integral

$(TARGET): $(OBJS)
	$(CXX) ./*.o $(INC_PATH) $(LIBS) $(CXXFLAGS) -o $(TARGET)
#	$(CXX) $(CXXFLAGS) ./*.o $(INC_PATH) -o $(TARGET)

$(OBJS): $(SRCFILES) 
	$(CXX) $(CXXFLAGS) -c $(SRCFILES) $(INC_PATH) $(LIBS)

#$(OBJS): $(SRCFILES) $(INC_PATH)/*.h
#	$(CXX) $<  $(CXXFLAGS) $(LIBS) -c -o $(TARGET)

clean:
	rm -rf *.out *.o $(OBJS) $(TARGET) 

