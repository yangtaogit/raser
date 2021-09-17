check: 
	python setup.py sdist 
	twine check dist/* 
	
upload: 
	twine upload dist/* 

clean: 
	rm -rf dist raser.egg-info  
 
