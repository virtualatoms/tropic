// document.addEventListener('DOMContentLoaded', function () {
window.addEventListener("load", function () {
    // Highlight active section based on scroll position
    window.addEventListener('scroll', function () {
        const headings = document.querySelectorAll('h1, h2, h3');
        const tocItems = document.querySelectorAll('.toc-item');

        let currentSection = '';
        headings.forEach(heading => {
            const rect = heading.getBoundingClientRect();
            if (rect.top <= 100) {
                currentSection = heading.id;
            }
        });

        tocItems.forEach(item => {
            if (item.getAttribute('href') === `#${currentSection}`) {
                item.classList.add('active');
            } else {
                item.classList.remove('active');
            }
        });
    });
});