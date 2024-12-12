# 使用的编译器
CC = mpic++

# 编译选项
CFLAGS = -std=c++11 -Iinclude

# 链接选项
LDFLAGS =

# 源文件和目标文件
SRCS = main.cpp lib/tmcmc.cpp
OBJS = $(SRCS:.cpp=.o)

# 可执行文件的名称
TARGET = tmcmc

# 默认目标，编译并生成可执行文件
all: $(TARGET)

# 链接目标文件生成最终可执行文件
$(TARGET): $(OBJS)
	@echo "Linking objects to create executable: $@"
	$(CC) $(LDFLAGS) -o $@ $^

# 编译源文件生成目标文件
%.o: %.cpp
	@echo "Compiling source file: $<"
	$(CC) $(CFLAGS) -c $< -o $@

# 清理生成的文件
clean:
	@echo "Cleaning up..."
	rm -f $(OBJS) $(TARGET)

# 打印调试信息
debug:
	@echo "Compiler: $(CC)"
	@echo "CFLAGS: $(CFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
	@echo "Source Files: $(SRCS)"
	@echo "Object Files: $(OBJS)"
	@echo "Target: $(TARGET)"

