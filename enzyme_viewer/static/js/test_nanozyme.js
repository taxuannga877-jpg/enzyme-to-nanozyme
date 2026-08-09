(function () {
    function testAPI() {
        const resultDiv = document.getElementById('result');
        resultDiv.textContent = 'Testing...';

        fetch('/api/list_nanozyme_types')
            .then(response => {
                console.log('Response status:', response.status);
                return response.json();
            })
            .then(data => {
                console.log('Response data:', data);
                resultDiv.textContent = '';
                const pre = document.createElement('pre');
                pre.textContent = JSON.stringify(data, null, 2);
                resultDiv.appendChild(pre);

                if (data.nanozyme_types && data.nanozyme_types.length > 0) {
                    const list = document.createElement('ul');
                    data.nanozyme_types.forEach(type => {
                        const li = document.createElement('li');
                        li.textContent = type;
                        list.appendChild(li);
                    });
                    resultDiv.appendChild(list);
                }
            })
            .catch(error => {
                console.error('Error:', error);
                resultDiv.textContent = 'Error: ' + error.message;
            });
    }

    window.testAPI = testAPI;
    window.addEventListener('load', () => {
        const button = document.getElementById('test-api-button');
        if (button) {
            button.addEventListener('click', testAPI);
        }
        testAPI();
    });
}());
