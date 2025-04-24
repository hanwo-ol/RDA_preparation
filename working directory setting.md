# working directory setting
가장 많이 오류가 나는 부분은 작업 디렉토리 설정입니다.(아마도요)

working directory를 위한 설명입니다.
* 저는 압축 파일을 바탕화면에 압축 해제를 한 상태입니다.

  
1. 우선, 실습 파일의 폴더가 있는 위치에 가시고,  
<img width="202" alt="작업 공간 세팅 1" src="https://github.com/user-attachments/assets/95d41d93-f9e1-4064-a9fb-6864d15521c7" />    

2. 우클릭 > 속성 클릭 하시면 됩니다.   
<img width="307" alt="image" src="https://github.com/user-attachments/assets/5f63a2b7-5d12-4ebf-be1b-c37437cef73f" />    

3. 이후에 해당 주소를 복사하시면 됩니다.   
<img width="303" alt="image" src="https://github.com/user-attachments/assets/2b2bade5-db6b-4ec5-8259-6530252f8abc" />

4. R studio로 돌아가셔서 붙여넣기 하시면 되는데...      
<img width="274" alt="image" src="https://github.com/user-attachments/assets/d035053b-daf1-42c9-a231-c198f736abf0" />

* 이대로 하면 아래처럼 오류가 뜹니다.
```
Error: '\U' used without hex digits in character string (<input>:1:15)
```
5. `\` 기호를 `/`로 바꿔주세요
6. 이후에... 



